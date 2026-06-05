# prepare_unmatched_ms2_domainmasst_merge.py
#
# Purpose:
#   Prepare MetOrigin-unmatched MS2 spectra for domainMASST/MASST-family searches,
#   optionally run the microbe_masst jobs.py-style workflow, parse JSON/TSV outputs,
#   and merge MASST evidence back into the MAPS + MetOriginDB table.
#
# Workflow:
#   MAPS_with_MetOriginDB_origin.csv + original MGF
#   -> remove metorigin_matched == TRUE feature IDs
#   -> extract unmatched feature MS2 spectra from MGF
#   -> rewrite TITLE as feature_id=<ID>|original_title=<old title>
#   -> write MetOrigin_unmatched_features_with_MS2.mgf
#   -> optionally run generated jobs_AHGMD_domainMASST.py
#   -> parse JSON/TSV outputs
#   -> merge evidence back to original MetOrigin table

import argparse
import json
import logging
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

YES_VALUES = {"yes", "y", "true", "t", "1"}


def clean_col(x: str) -> str:
    x = x.strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def is_yes(x) -> bool:
    if pd.isna(x):
        return False
    return str(x).strip().lower() in YES_VALUES


def extract_feature_id_from_text(text: str) -> Optional[str]:
    if not text:
        return None
    patterns = [
        r"feature_id[=:| _-]*([A-Za-z0-9_.-]+)",
        r"feature\.id[=:| _-]*([A-Za-z0-9_.-]+)",
        r"featureID[=:| _-]*([A-Za-z0-9_.-]+)",
        r"Feature[_ -]?([0-9]+)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return str(m.group(1))
    return None


def parse_mgf_blocks(mgf_path: Path) -> List[Dict]:
    lines = mgf_path.read_text(errors="replace").splitlines()
    blocks = []
    current = []
    for line in lines:
        if line.strip().upper() == "BEGIN IONS":
            current = [line]
        elif line.strip().upper() == "END IONS":
            current.append(line)
            blocks.append(parse_one_mgf_block(current))
            current = []
        elif current:
            current.append(line)
    return blocks


def get_mgf_field(block_lines: List[str], key: str) -> Optional[str]:
    key_upper = key.upper() + "="
    for line in block_lines:
        if line.upper().startswith(key_upper):
            return line.split("=", 1)[1].strip()
    return None


def set_mgf_field(block_lines: List[str], key: str, value: str) -> List[str]:
    key_upper = key.upper() + "="
    out = []
    replaced = False
    for line in block_lines:
        if line.upper().startswith(key_upper):
            out.append(f"{key}={value}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        insert_at = 1 if out and out[0].strip().upper() == "BEGIN IONS" else 0
        out.insert(insert_at, f"{key}={value}")
    return out


def parse_one_mgf_block(block_lines: List[str]) -> Dict:
    title = get_mgf_field(block_lines, "TITLE")
    scans = get_mgf_field(block_lines, "SCANS")
    feature_id = extract_feature_id_from_text(title or "") or extract_feature_id_from_text(scans or "") or scans
    return {
        "feature_id": str(feature_id) if feature_id is not None else None,
        "title": title,
        "scans": scans,
        "pepmass": get_mgf_field(block_lines, "PEPMASS"),
        "rtinseconds": get_mgf_field(block_lines, "RTINSECONDS"),
        "charge": get_mgf_field(block_lines, "CHARGE"),
        "block_lines": block_lines,
    }


def prepare_unmatched_mgf(
    metorigin_csv: Path,
    input_mgf: Path,
    out_dir: Path,
    metorigin_matched_col: str = "metorigin_matched",
) -> Tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(metorigin_csv)
    df.columns = [clean_col(c) for c in df.columns]
    if "feature_id" not in df.columns:
        raise ValueError("MetOrigin table must contain feature_id.")
    if metorigin_matched_col not in df.columns:
        raise ValueError(f"MetOrigin table must contain {metorigin_matched_col}.")

    df["feature_id"] = df["feature_id"].astype(str)
    matched_ids = set(df.loc[df[metorigin_matched_col].apply(is_yes), "feature_id"].astype(str))
    unmatched_df = df.loc[~df["feature_id"].isin(matched_ids)].copy()
    unmatched_ids = set(unmatched_df["feature_id"].astype(str))

    logger.info("Total MetOrigin rows: %d", len(df))
    logger.info("MetOrigin-matched feature IDs removed: %d", len(matched_ids))
    logger.info("MetOrigin-unmatched feature IDs retained: %d", len(unmatched_ids))

    mgf_blocks = parse_mgf_blocks(input_mgf)
    kept_blocks = []
    manifest_rows = []

    for b in mgf_blocks:
        fid = b["feature_id"]
        if fid is None:
            continue
        fid = str(fid)
        if fid not in unmatched_ids:
            continue
        original_title = b["title"] or ""
        new_title = f"feature_id={fid}|original_title={original_title}"
        new_lines = set_mgf_field(b["block_lines"], "TITLE", new_title)
        new_lines = set_mgf_field(new_lines, "SCANS", fid)
        kept_blocks.append(new_lines)
        manifest_rows.append({
            "feature_id": fid,
            "original_title": original_title,
            "new_title": new_title,
            "scans": fid,
            "pepmass": b["pepmass"],
            "rtinseconds": b["rtinseconds"],
            "charge": b["charge"],
        })

    manifest = pd.DataFrame(manifest_rows).drop_duplicates()
    if manifest.empty:
        logger.warning("No unmatched MS2 spectra were found in the MGF. Check feature_id format in MGF TITLE/SCANS.")
        unmatched_with_ms2 = unmatched_df.iloc[0:0].copy()
    else:
        unmatched_with_ms2 = unmatched_df.merge(
            manifest.groupby("feature_id").agg(
                n_ms2_spectra=("feature_id", "size"),
                masst_query_titles=("new_title", lambda x: " || ".join(pd.unique(x.astype(str)))),
                precursor_mz_mgf=("pepmass", lambda x: " || ".join(pd.unique(x.dropna().astype(str)))),
                rt_seconds_mgf=("rtinseconds", lambda x: " || ".join(pd.unique(x.dropna().astype(str)))),
            ).reset_index(),
            on="feature_id",
            how="inner",
        )

    out_mgf = out_dir / "MetOrigin_unmatched_features_with_MS2.mgf"
    out_manifest = out_dir / "MetOrigin_unmatched_features_with_MS2_manifest.csv"
    out_unmatched_table = out_dir / "MetOrigin_unmatched_features_with_MS2.csv"

    with out_mgf.open("w", newline="\n") as f:
        for block in kept_blocks:
            f.write("\n".join(block))
            f.write("\n")
    manifest.to_csv(out_manifest, index=False)
    unmatched_with_ms2.to_csv(out_unmatched_table, index=False)

    logger.info("Unmatched MS2 spectra exported: %d", len(kept_blocks))
    logger.info("Unmatched features with MS2: %d", unmatched_with_ms2["feature_id"].nunique() if "feature_id" in unmatched_with_ms2 else 0)
    return out_mgf, out_manifest, out_unmatched_table


def write_domainmasst_job_py(
    out_dir: Path,
    mgf_file: Path,
    output_prefix: Path,
    analog: bool = False,
    min_cos: float = 0.7,
    mz_tol: float = 0.02,
    precursor_mz_tol: float = 0.02,
    min_matched_signals: int = 3,
    parallel_queries: int = 5,
) -> Path:
    job_py = out_dir / "jobs_AHGMD_domainMASST.py"
    project_dir = Path.cwd().resolve()
    upstream_code_dir = project_dir / "external" / "microbe_masst" / "code"
    script = f'''import logging
import sys
from pathlib import Path

UPSTREAM_CODE_DIR = Path(r"{upstream_code_dir}")
if UPSTREAM_CODE_DIR.exists():
    sys.path.insert(0, str(UPSTREAM_CODE_DIR))

import masst_batch_client
import masst_utils

logging.basicConfig(
    level=logging.DEBUG,
    handlers=[
        logging.FileHandler(Path(r"{(out_dir / 'domainmasst_batch.log').resolve()}"), mode="w", encoding="utf-8"),
        logging.StreamHandler(),
    ],
    force=True,
)
logger = logging.getLogger(__name__)

files = [
    (r"{mgf_file}", r"{output_prefix}")
]

if __name__ == "__main__":
    failures = 0
    for file, out_file in files:
        try:
            logger.info("Starting new job for input: {{}}".format(file))
            sep = "," if file.endswith("csv") else "\\t"
            success_rate = masst_batch_client.run_on_usi_list_or_mgf_file(
                in_file=file,
                out_file_no_extension=out_file,
                min_cos={min_cos},
                mz_tol={mz_tol},
                precursor_mz_tol={precursor_mz_tol},
                min_matched_signals={min_matched_signals},
                database=masst_utils.DataBase.metabolomicspanrepo_index_nightly,
                parallel_queries={parallel_queries},
                skip_existing=True,
                analog={str(analog)},
                sep=sep,
            )
            logger.info("domainMASST success rate: %.3f", success_rate)
            if success_rate < 1:
                failures += 1
        except Exception as e:
            logger.exception(e)
            failures += 1
    sys.exit(1 if failures else 0)
'''
    job_py.write_text(script)
    logger.info("Wrote domainMASST job script: %s", job_py)
    return job_py


def run_job_py(job_py: Path, cwd: Optional[Path] = None) -> None:
    logger.info("Running: python %s", job_py)
    subprocess.run([sys.executable, str(job_py.resolve())], cwd=str(cwd) if cwd else None, check=True)


def flatten_json_summary(obj, prefix="") -> Dict[str, str]:
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            new_prefix = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, (dict, list)):
                out.update(flatten_json_summary(v, new_prefix))
            else:
                out[new_prefix] = v
    elif isinstance(obj, list):
        out[f"{prefix}.list_length"] = len(obj)
        text_hits = []
        for item in obj[:20]:
            if isinstance(item, dict):
                label = item.get("name") or item.get("label") or item.get("id") or item.get("taxon") or item.get("ontology")
                if label:
                    text_hits.append(str(label))
            else:
                text_hits.append(str(item))
        if text_hits:
            out[f"{prefix}.top_items"] = " || ".join(text_hits)
    return out


def parse_feature_id_from_json_or_filename(json_path: Path, data) -> Optional[str]:
    text_candidates = [json_path.name, json_path.stem]
    def collect_strings(x, max_n=100):
        vals = []
        if isinstance(x, dict):
            for k, v in x.items():
                if k.lower() in {"title", "query", "query_title", "usi", "spectrum", "scan", "scans", "id", "name"}:
                    vals.append(str(v))
                vals.extend(collect_strings(v, max_n=max_n))
        elif isinstance(x, list):
            for item in x[:20]:
                vals.extend(collect_strings(item, max_n=max_n))
        elif isinstance(x, str):
            vals.append(x)
        return vals[:max_n]
    text_candidates.extend(collect_strings(data))
    for txt in text_candidates:
        fid = extract_feature_id_from_text(txt)
        if fid:
            return str(fid)
    return None


def load_manifest_feature_ids(manifest_csv: Path) -> List[str]:
    if not manifest_csv.exists():
        return []
    manifest = pd.read_csv(manifest_csv)
    if manifest.empty or "feature_id" not in manifest.columns:
        return []
    return sorted(manifest["feature_id"].dropna().astype(str).unique(), key=len, reverse=True)


def feature_id_from_output_path(path: Path, feature_ids: List[str]) -> Optional[str]:
    fid = extract_feature_id_from_text(path.name)
    if fid:
        return str(fid)
    stem = path.stem
    for suffix in [
        "_unfiltered_matches",
        "_analog_matches",
        "_count_domain",
        "_matches",
        "_library",
        "_datasets",
        "_microbe",
        "_plant",
        "_tissue",
        "_food",
        "_microbiome",
        "_combined",
        "_domain",
    ]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for feature_id in feature_ids:
        if re.search(rf"(^|_){re.escape(feature_id)}$", stem):
            return feature_id
    return None


def parse_domainmasst_outputs(domainmasst_dir: Path, manifest_csv: Path) -> pd.DataFrame:
    manifest_feature_ids = load_manifest_feature_ids(manifest_csv)
    rows = []
    for json_path in domainmasst_dir.rglob("*.json"):
        try:
            data = json.loads(json_path.read_text(errors="replace"))
        except Exception as e:
            logger.warning("Could not parse JSON %s: %s", json_path, e)
            continue
        fid = parse_feature_id_from_json_or_filename(json_path, data) or feature_id_from_output_path(json_path, manifest_feature_ids)
        domain = "unknown"
        lower_name = json_path.name.lower()
        for d in ["microbe", "plant", "tissue", "food", "microbiome", "combined"]:
            if d in lower_name:
                domain = d
                break
        flat = flatten_json_summary(data)
        rows.append({
            "feature_id": fid,
            "masst_json_file": str(json_path),
            "masst_domain": domain,
            "json_top_items": " || ".join(str(v) for k, v in flat.items() if k.endswith("top_items") and pd.notna(v))[:2000],
            "json_flat_key_count": len(flat),
        })
    json_summary = pd.DataFrame(rows)
    if len(json_summary) == 0:
        logger.warning("No JSON files found or parsed in %s", domainmasst_dir)
        return pd.DataFrame(columns=["feature_id", "masst_json_files", "masst_domains_detected", "masst_json_top_items"])
    if json_summary["feature_id"].isna().any():
        logger.warning("Some JSON files did not contain feature_id. They cannot be merged unless filenames/query titles preserve feature_id.")
    json_summary = json_summary.dropna(subset=["feature_id"]).copy()
    json_summary["feature_id"] = json_summary["feature_id"].astype(str)
    return json_summary.groupby("feature_id").agg(
        masst_json_files=("masst_json_file", lambda x: " || ".join(pd.unique(x.astype(str)))),
        masst_domains_detected=("masst_domain", lambda x: " || ".join(pd.unique(x.astype(str)))),
        masst_json_top_items=("json_top_items", lambda x: " || ".join([i for i in pd.unique(x.astype(str)) if i and i != "nan"])[:3000]),
    ).reset_index()


def parse_tsv_outputs(domainmasst_dir: Path, manifest_csv: Path) -> pd.DataFrame:
    manifest_feature_ids = load_manifest_feature_ids(manifest_csv)
    rows = []
    for tsv_path in domainmasst_dir.rglob("*.tsv"):
        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except Exception:
            continue
        file_type = "tsv"
        name = tsv_path.name.lower()
        if "_matches" in name:
            file_type = "matches"
        elif "_library" in name:
            file_type = "library"
        elif "_datasets" in name:
            file_type = "datasets"
        elif "_count_domain" in name:
            file_type = "count_domain"
        fid = feature_id_from_output_path(tsv_path, manifest_feature_ids)
        if fid is None:
            for col in df.columns:
                if col.lower() in {"title", "query", "query_title", "scan", "scans", "spectrum_id", "usi"}:
                    candidate = df[col].dropna().astype(str).map(extract_feature_id_from_text).dropna()
                    if len(candidate):
                        fid = str(candidate.iloc[0])
                        break
        rows.append({
            "feature_id": fid,
            "masst_tsv_file": str(tsv_path),
            "masst_tsv_type": file_type,
            "masst_tsv_n_rows": len(df),
            "masst_tsv_columns": ";".join(df.columns.astype(str)),
        })
    tsv_summary = pd.DataFrame(rows)
    if len(tsv_summary) == 0:
        return pd.DataFrame(columns=["feature_id", "masst_tsv_files", "masst_n_matches_rows", "masst_n_library_rows", "masst_n_dataset_rows", "masst_n_count_domain_rows"])
    tsv_summary = tsv_summary.dropna(subset=["feature_id"]).copy()
    tsv_summary["feature_id"] = tsv_summary["feature_id"].astype(str)
    pivot = tsv_summary.pivot_table(index="feature_id", columns="masst_tsv_type", values="masst_tsv_n_rows", aggfunc="sum", fill_value=0).reset_index()
    for col in ["matches", "library", "datasets", "count_domain"]:
        if col not in pivot.columns:
            pivot[col] = 0
    files = tsv_summary.groupby("feature_id").agg(masst_tsv_files=("masst_tsv_file", lambda x: " || ".join(pd.unique(x.astype(str))))).reset_index()
    out = files.merge(pivot, on="feature_id", how="left")
    return out.rename(columns={"matches": "masst_n_matches_rows", "library": "masst_n_library_rows", "datasets": "masst_n_dataset_rows", "count_domain": "masst_n_count_domain_rows"})


def merge_back_to_metorigin(metorigin_csv: Path, manifest_csv: Path, domainmasst_dir: Path, out_csv: Path) -> Path:
    metorigin = pd.read_csv(metorigin_csv)
    metorigin.columns = [clean_col(c) for c in metorigin.columns]
    if "feature_id" not in metorigin.columns:
        raise ValueError("MetOrigin table must contain feature_id.")
    metorigin["feature_id"] = metorigin["feature_id"].astype(str)
    manifest = pd.read_csv(manifest_csv)
    if manifest.empty:
        searched = pd.DataFrame(columns=["feature_id", "masst_searched", "masst_n_submitted_spectra", "masst_query_titles"])
    else:
        manifest["feature_id"] = manifest["feature_id"].astype(str)
    searched = manifest.groupby("feature_id").agg(
            masst_searched=("feature_id", lambda x: True),
            masst_n_submitted_spectra=("feature_id", "size"),
            masst_query_titles=("new_title", lambda x: " || ".join(pd.unique(x.astype(str)))),
    ).reset_index()
    json_summary = parse_domainmasst_outputs(domainmasst_dir, manifest_csv)
    tsv_summary = parse_tsv_outputs(domainmasst_dir, manifest_csv)
    merged = metorigin.merge(searched, on="feature_id", how="left").merge(json_summary, on="feature_id", how="left").merge(tsv_summary, on="feature_id", how="left")
    merged["masst_searched"] = merged["masst_searched"].fillna(False)
    merged["origin_evidence_merge_note"] = merged["masst_searched"].map(
        lambda x: "MetOrigin identity-based evidence plus domainMASST spectrum/repository evidence kept as separate evidence blocks." if x else "No domainMASST search for this feature, either MetOrigin-matched already or no MS2 spectrum was exported."
    )
    merged.to_csv(out_csv, index=False)
    logger.info("Merged output written: %s", out_csv)
    return out_csv


def main():
    parser = argparse.ArgumentParser(description="Prepare MetOrigin-unmatched MS2 MGF, run domainMASST jobs.py-style search, and merge JSON/TSV evidence back.")
    parser.add_argument("--metorigin_csv", required=True, help="MAPS + MetOriginDB result CSV.")
    parser.add_argument("--mgf", required=True, help="Original MGF containing feature_id in TITLE or SCANS.")
    parser.add_argument("--out_dir", default="Graph1_MASST_unmatched_MS2_outputs")
    parser.add_argument("--run_masst", action="store_true", help="Actually run the generated domainMASST job script.")
    parser.add_argument("--domainmasst_dir", default=None, help="Existing domainMASST output directory if already run.")
    parser.add_argument("--analog", action="store_true", help="Run analog MASST mode. Default is exact mode.")
    parser.add_argument("--min_cos", type=float, default=0.7)
    parser.add_argument("--mz_tol", type=float, default=0.02)
    parser.add_argument("--precursor_mz_tol", type=float, default=0.02)
    parser.add_argument("--min_matched_signals", type=int, default=3)
    parser.add_argument("--parallel_queries", type=int, default=5)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    domainmasst_dir = Path(args.domainmasst_dir) if args.domainmasst_dir else out_dir / "domainMASST"
    domainmasst_dir.mkdir(parents=True, exist_ok=True)

    out_mgf, manifest_csv, unmatched_table = prepare_unmatched_mgf(Path(args.metorigin_csv), Path(args.mgf), out_dir)
    output_prefix = domainmasst_dir / "AHGMD_unmatched"
    job_py = write_domainmasst_job_py(out_dir, out_mgf.resolve(), output_prefix.resolve(), analog=args.analog, min_cos=args.min_cos, mz_tol=args.mz_tol, precursor_mz_tol=args.precursor_mz_tol, min_matched_signals=args.min_matched_signals, parallel_queries=args.parallel_queries)
    if args.run_masst:
        upstream_code_dir = Path("external/microbe_masst/code").resolve()
        run_job_py(job_py, cwd=upstream_code_dir if upstream_code_dir.exists() else None)
    merged_csv = out_dir / "MAPS_with_MetOriginDB_plus_domainMASST_origin_evidence.csv"
    merge_back_to_metorigin(Path(args.metorigin_csv), manifest_csv, domainmasst_dir, merged_csv)

    logger.info("Done.")
    logger.info("Unmatched MGF: %s", out_mgf)
    logger.info("Manifest: %s", manifest_csv)
    logger.info("Unmatched table: %s", unmatched_table)
    logger.info("Job script: %s", job_py)
    logger.info("Merged table: %s", merged_csv)


if __name__ == "__main__":
    main()
