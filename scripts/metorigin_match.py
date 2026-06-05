#!/usr/bin/env python
"""
MAPS -> MetOriginDB matching only.

Purpose
-------
Join MAPS final-annotation-df.csv to MetOriginDB.csv while preserving one row per
MAPS feature/annotation.

Important MetOriginDB behavior
------------------------------
MetOriginDB can contain multiple rows for the same metabolite, PubChem CID, SMILES,
HMDB ID, etc. This script aggregates all MetOriginDB rows per matching key before
joining to MAPS.

Main mappings
-------------
MAPS CID / cid       -> MetOriginDB PUBCHEM_COMPOUND_ID
MAPS smiles / SMILES -> MetOriginDB SMILES_ID

Matching priority
-----------------
1. PubChem CID
2. SMILES exact
3. InChIKey
4. HMDB ID
5. KEGG ID
6. compound name exact

Run
---
python scripts/metorigin_match.py --config config/config.yaml
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import pandas as pd
import yaml


def setup_logging(log_file: Optional[Path] = None) -> None:
    handlers = [logging.StreamHandler()]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def clean_colname(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def clean_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cleaned_cols = [clean_colname(c) for c in df.columns]
    seen: Dict[str, int] = {}
    unique_cols = []
    for c in cleaned_cols:
        if c not in seen:
            seen[c] = 0
            unique_cols.append(c)
        else:
            seen[c] += 1
            unique_cols.append(f"{c}_{seen[c]}")
    df.columns = unique_cols
    return df


def clean_id_series(s: pd.Series) -> pd.Series:
    out = s.astype("string").str.strip()
    out = out.str.replace(r"\.0$", "", regex=True)
    out = out.replace({
        "": pd.NA, "NA": pd.NA, "N/A": pd.NA, "na": pd.NA, "n/a": pd.NA,
        "nan": pd.NA, "NaN": pd.NA, "None": pd.NA, "none": pd.NA,
        "NULL": pd.NA, "null": pd.NA,
    })
    return out


def clean_name_series(s: pd.Series, case_insensitive: bool = True) -> pd.Series:
    out = s.astype("string").str.replace(r"\s+", " ", regex=True).str.strip()
    out = out.replace({
        "": pd.NA, "NA": pd.NA, "N/A": pd.NA, "na": pd.NA, "n/a": pd.NA,
        "nan": pd.NA, "NaN": pd.NA, "None": pd.NA, "none": pd.NA,
        "NULL": pd.NA, "null": pd.NA,
    })
    if case_insensitive:
        out = out.str.lower()
    return out


def pick_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    for c in [clean_colname(x) for x in candidates]:
        if c in df.columns:
            return c
    return None


def apply_overrides(cols: Dict[str, Optional[str]], overrides: Dict[str, Any], df: pd.DataFrame, label: str) -> Dict[str, Optional[str]]:
    if not overrides:
        return cols
    for logical_name, original_col in overrides.items():
        logical_name_clean = clean_colname(logical_name)
        original_col_clean = clean_colname(str(original_col))
        if original_col_clean in df.columns:
            cols[logical_name_clean] = original_col_clean
            logging.info("%s override: %s -> %s", label, logical_name_clean, original_col_clean)
        else:
            logging.warning("%s override requested %s -> %s, but column was not found.", label, logical_name_clean, original_col_clean)
    return cols


def get_nested(config: Dict[str, Any], possible_paths: Iterable[Tuple[str, ...]], default: Any = None) -> Any:
    for path in possible_paths:
        current: Any = config
        ok = True
        for key in path:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                ok = False
                break
        if ok:
            return current
    return default


def resolve_path(path_value: str, base_dir: Path) -> Path:
    p = Path(path_value)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def detect_maps_columns(maps: pd.DataFrame) -> Dict[str, Optional[str]]:
    cols = {
        "feature_id": pick_col(maps, ["feature_id", "featureid", "id"]),
        "usi": pick_col(maps, ["feature_usi", "usi"]),
        "name": pick_col(maps, ["compound_name", "compound.name", "name", "title"]),
        "smiles": pick_col(maps, ["smiles", "SMILES", "canonical_smiles"]),
        "inchikey": pick_col(maps, ["inchikey", "inchi_key", "inchi.key"]),
        "cid": pick_col(maps, ["CID", "cid", "pubchem_cid", "pubchem_id", "pubchem_compound_id"]),
        "hmdb": pick_col(maps, ["hmdb_id", "primary_hmdb_id", "hmdb"]),
        "kegg": pick_col(maps, ["kegg_id", "kegg"]),
        "confidence": pick_col(maps, ["confidence_level", "confidence.level", "annotation_confidence_level"]),
        "annotation_type": pick_col(maps, ["annotation_type", "annotation.type"]),
        "mz": pick_col(maps, ["mz"]),
        "rt": pick_col(maps, ["rt"]),
        "formula": pick_col(maps, ["formula"]),
        "iupac": pick_col(maps, ["iupac"]),
        "mono_mass": pick_col(maps, ["monoisotopic_mass", "monoisotopic.mass", "mass"]),
    }
    if "cid" in maps.columns:
        cols["cid"] = "cid"
    if "smiles" in maps.columns:
        cols["smiles"] = "smiles"
    return cols


def detect_metorigin_columns(metorigin: pd.DataFrame) -> Dict[str, Optional[str]]:
    cols = {
        "name": pick_col(metorigin, ["compound_name", "compound_name_id", "compound.name", "name", "metabolite_name", "title", "chemical_name"]),
        "smiles": pick_col(metorigin, ["SMILES_ID", "smiles_id", "smiles", "canonical_smiles", "isomeric_smiles", "pubchem_smiles"]),
        "inchikey": pick_col(metorigin, ["inchikey", "inchi_key", "inchi.key", "standard_inchikey"]),
        "cid": pick_col(metorigin, ["PUBCHEM_COMPOUND_ID", "pubchem_compound_id", "cid", "pubchem_cid", "pubchem_id", "pubchem"]),
        "hmdb": pick_col(metorigin, ["hmdb_id", "primary_hmdb_id", "hmdb", "hmdbid"]),
        "kegg": pick_col(metorigin, ["kegg_id", "kegg", "keggid"]),
    }
    if "pubchem_compound_id" in metorigin.columns:
        cols["cid"] = "pubchem_compound_id"
    if "smiles_id" in metorigin.columns:
        cols["smiles"] = "smiles_id"
    return cols


ORIGIN_FLAG_COLUMNS = [
    "from_human", "from_bacteria", "from_plant", "from_animal",
    "from_environment", "from_drug", "from_food",
]

ORIGIN_DETAIL_COLUMNS = [
    "from_which_part", "from_which_bacteria", "bacteria_ncbi_id", "bacteria_phylum",
    "bacteria_class", "bacteria_order", "bacteria_family", "bacteria_genus",
    "bacteria_species", "from_which_plant", "from_which_animal",
    "from_which_environment", "from_which_drug", "from_which_food",
]

ORIGIN_COLUMNS_TO_KEEP = ORIGIN_FLAG_COLUMNS + ORIGIN_DETAIL_COLUMNS

METORIGIN_IDENTITY_COLUMNS = [
    "pubchem_compound_id", "smiles_id", "compound_name", "compound_name_id",
    "hmdb_id", "kegg_id", "inchikey", "mimedb_id", "lotus_id",
]


def add_match_key(df: pd.DataFrame, source_col: Optional[str], key_col: str, kind: str, case_insensitive_names: bool = True) -> pd.DataFrame:
    df = df.copy()
    if source_col is None or source_col not in df.columns:
        return df
    if kind == "name":
        df[key_col] = clean_name_series(df[source_col], case_insensitive=case_insensitive_names)
    else:
        df[key_col] = clean_id_series(df[source_col])
    return df


def prepare_match_keys(maps: pd.DataFrame, metorigin: pd.DataFrame, maps_cols: Dict[str, Optional[str]], met_cols: Dict[str, Optional[str]], case_insensitive_names: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    maps = maps.copy()
    metorigin = metorigin.copy()
    key_specs = [
        ("cid", "match_cid", "id"),
        ("smiles", "match_smiles", "id"),
        ("inchikey", "match_inchikey", "id"),
        ("hmdb", "match_hmdb", "id"),
        ("kegg", "match_kegg", "id"),
        ("name", "match_name", "name"),
    ]
    for logical, key_col, kind in key_specs:
        maps = add_match_key(maps, maps_cols.get(logical), key_col, kind, case_insensitive_names=case_insensitive_names)
        metorigin = add_match_key(metorigin, met_cols.get(logical), key_col, kind, case_insensitive_names=case_insensitive_names)
    return maps, metorigin


def collapse_yes_unknown(values: pd.Series) -> str:
    vals = (
        values.dropna().astype(str).str.strip().replace("", pd.NA).dropna().str.lower().tolist()
    )
    if any(v == "yes" for v in vals):
        return "Yes"
    if any(v == "unknown" for v in vals):
        return "Unknown"
    if any(v == "no" for v in vals):
        return "No"
    return ""


def collapse_unique_text(values: pd.Series, sep: str = " | ") -> str:
    vals = values.dropna().astype(str).str.strip().replace("", pd.NA).dropna()
    if vals.empty:
        return ""
    unique_vals = sorted(set(vals.tolist()))
    return sep.join(unique_vals)


def build_metorigin_subset(met_db: pd.DataFrame, key: str) -> pd.DataFrame:
    """
    Build a one-row-per-key MetOriginDB lookup table.

    MetOriginDB can contain multiple rows for the same key, so this aggregates
    all rows rather than keeping the first row.
    """
    keep_cols = [key]
    for c in ORIGIN_COLUMNS_TO_KEEP:
        if c in met_db.columns and c not in keep_cols:
            keep_cols.append(c)
    for c in METORIGIN_IDENTITY_COLUMNS:
        if c in met_db.columns and c not in keep_cols:
            keep_cols.append(c)

    met_tmp = met_db.loc[:, keep_cols].dropna(subset=[key]).copy()
    if met_tmp.empty:
        return met_tmp

    agg_dict: Dict[str, Any] = {}
    for c in ORIGIN_FLAG_COLUMNS:
        if c in met_tmp.columns:
            agg_dict[c] = collapse_yes_unknown
    for c in ORIGIN_DETAIL_COLUMNS:
        if c in met_tmp.columns:
            agg_dict[c] = collapse_unique_text
    for c in METORIGIN_IDENTITY_COLUMNS:
        if c in met_tmp.columns and c != key:
            agg_dict[c] = collapse_unique_text

    met_small = met_tmp.groupby(key, dropna=False).agg(agg_dict).reset_index()
    record_counts = met_tmp.groupby(key, dropna=False).size().reset_index(name="metorigin_n_database_records")
    met_small = met_small.merge(record_counts, on=key, how="left")

    rename_map = {c: f"metorigin_{c}" for c in METORIGIN_IDENTITY_COLUMNS if c in met_small.columns and c != key}
    met_small = met_small.rename(columns=rename_map)
    return met_small


def join_by_key(unmatched: pd.DataFrame, met_db: pd.DataFrame, key: str, match_method: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if key not in unmatched.columns or key not in met_db.columns:
        logging.info("Skipping %s because key column is missing.", match_method)
        return pd.DataFrame(), unmatched

    u = unmatched[unmatched[key].notna()].copy()
    no_key = unmatched[unmatched[key].isna()].copy()
    if u.empty:
        logging.info("Skipping %s because no MAPS rows have this key.", match_method)
        return pd.DataFrame(), unmatched

    met_small = build_metorigin_subset(met_db, key)
    if met_small.empty:
        logging.info("Skipping %s because no MetOriginDB rows have this key.", match_method)
        return pd.DataFrame(), unmatched

    joined = u.merge(met_small, how="left", on=key, indicator=True)
    matched = joined[joined["_merge"] == "both"].copy().drop(columns=["_merge"])
    still_unmatched = joined[joined["_merge"] == "left_only"].copy().drop(columns=["_merge"])

    if not matched.empty:
        matched["metorigin_matched"] = True
        matched["metorigin_match_method"] = match_method

    still_unmatched = still_unmatched.loc[:, [c for c in unmatched.columns if c in still_unmatched.columns]]
    next_unmatched = pd.concat([still_unmatched, no_key], ignore_index=True)
    logging.info("Matched by %s: %s", match_method, len(matched))
    return matched, next_unmatched


def run_matching(maps: pd.DataFrame, metorigin: pd.DataFrame, maps_cols: Dict[str, Optional[str]], met_cols: Dict[str, Optional[str]], case_insensitive_names: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    maps2, met2 = prepare_match_keys(maps, metorigin, maps_cols, met_cols, case_insensitive_names=case_insensitive_names)
    match_plan = [
        ("match_cid", "PubChem CID"),
        ("match_smiles", "SMILES exact"),
        ("match_inchikey", "InChIKey"),
        ("match_hmdb", "HMDB ID"),
        ("match_kegg", "KEGG ID"),
        ("match_name", "compound name exact"),
    ]
    unmatched = maps2.copy()
    matched_parts = []
    for key, method in match_plan:
        matched, unmatched = join_by_key(unmatched, met2, key, method)
        if not matched.empty:
            matched_parts.append(matched)

    unmatched = unmatched.copy()
    unmatched["metorigin_matched"] = False
    unmatched["metorigin_match_method"] = "unmatched"
    final = pd.concat(matched_parts + [unmatched], ignore_index=True) if matched_parts else unmatched
    summary = (
        final.groupby(["metorigin_matched", "metorigin_match_method"], dropna=False)
        .size()
        .reset_index(name="n_features")
        .sort_values(["metorigin_matched", "n_features"], ascending=[False, False])
    )
    return final, summary


def is_yes_value(x: object) -> bool:
    if pd.isna(x):
        return False
    return str(x).strip().lower() in {"yes", "y", "true", "t", "1"}


def add_origin_summary_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    category_map = {
        "from_human": "human",
        "from_bacteria": "bacteria",
        "from_plant": "plant",
        "from_animal": "animal",
        "from_environment": "environment",
        "from_drug": "drug",
        "from_food": "food",
    }
    for col in ORIGIN_FLAG_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA

    def summarize_row(row: pd.Series) -> str:
        cats = []
        for col, label in category_map.items():
            if is_yes_value(row.get(col, pd.NA)):
                cats.append(label)
        return ";".join(cats) if cats else ""

    df["metorigin_origin_categories"] = df.apply(summarize_row, axis=1)
    df["metorigin_n_origin_categories"] = df["metorigin_origin_categories"].apply(lambda x: 0 if x == "" else len(str(x).split(";")))
    return df


def drop_internal_match_columns(df: pd.DataFrame) -> pd.DataFrame:
    internal_cols = ["match_cid", "match_smiles", "match_inchikey", "match_hmdb", "match_kegg", "match_name"]
    return df.drop(columns=[c for c in internal_cols if c in df.columns], errors="ignore")


def make_compact_output(df: pd.DataFrame) -> pd.DataFrame:
    df = add_origin_summary_columns(df)
    df = drop_internal_match_columns(df)
    preferred_maps_cols = [
        "feature_id", "feature_usi", "usi", "mz", "rt", "compound_name", "smiles",
        "formula", "iupac", "monoisotopic_mass", "cid", "hmdb_id", "annotation_type",
        "confidence_level", "confidence_score",
    ]
    match_cols = [
        "metorigin_matched", "metorigin_match_method", "metorigin_n_database_records",
        "metorigin_origin_categories", "metorigin_n_origin_categories",
    ]
    metorigin_identity_cols = [
        "metorigin_pubchem_compound_id", "metorigin_smiles_id", "metorigin_compound_name",
        "metorigin_compound_name_id", "metorigin_hmdb_id", "metorigin_kegg_id",
        "metorigin_inchikey", "metorigin_mimedb_id", "metorigin_lotus_id",
    ]
    compact_cols = []
    for c in preferred_maps_cols + match_cols + ORIGIN_FLAG_COLUMNS + metorigin_identity_cols + ORIGIN_DETAIL_COLUMNS:
        if c in df.columns and c not in compact_cols:
            compact_cols.append(c)
    remaining_cols = [c for c in df.columns if c not in compact_cols and not c.startswith("match_")]
    return df.loc[:, compact_cols + remaining_cols]


def make_excel_friendly_preview(df: pd.DataFrame, max_detail_chars: int = 300) -> pd.DataFrame:
    out = df.copy()
    for c in ORIGIN_DETAIL_COLUMNS:
        if c in out.columns:
            out[c] = out[c].astype("string").apply(
                lambda x: x if pd.isna(x) or len(str(x)) <= max_detail_chars else str(x)[:max_detail_chars] + "..."
            )
    return out


def write_column_detection_report(out_file: Path, maps: pd.DataFrame, metorigin: pd.DataFrame, maps_cols: Dict[str, Optional[str]], met_cols: Dict[str, Optional[str]]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8") as f:
        f.write("MAPS -> MetOriginDB column detection report\n")
        f.write("=" * 80 + "\n\n")
        f.write("Detected MAPS logical columns\n")
        f.write("-" * 80 + "\n")
        for k, v in maps_cols.items():
            f.write(f"{k}: {v}\n")
        f.write("\nDetected MetOriginDB logical columns\n")
        f.write("-" * 80 + "\n")
        for k, v in met_cols.items():
            f.write(f"{k}: {v}\n")
        f.write("\nAll MAPS columns after cleaning\n")
        f.write("-" * 80 + "\n")
        for c in maps.columns:
            f.write(f"{c}\n")
        f.write("\nAll MetOriginDB columns after cleaning\n")
        f.write("-" * 80 + "\n")
        for c in metorigin.columns:
            f.write(f"{c}\n")


def load_config(config_path: Path) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    if config is None:
        config = {}
    if not isinstance(config, dict):
        raise ValueError("Config file must contain a YAML dictionary.")
    return config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config/config.yaml")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    project_dir = config_path.parent.parent
    config = load_config(config_path)

    maps_file = get_nested(config, [("maps_file",), ("input", "maps_file"), ("inputs", "maps_file"), ("paths", "maps_file"), ("files", "maps_file")], default="data/maps/final-annotation-df.csv")
    metorigin_file = get_nested(config, [("metorigin_file",), ("input", "metorigin_file"), ("inputs", "metorigin_file"), ("paths", "metorigin_file"), ("files", "metorigin_file")], default="data/metorigin/MetOriginDB.csv")
    output_file = get_nested(config, [("output_file",), ("outputs", "matched_output"), ("outputs", "output_file"), ("paths", "output_file")], default="outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv")
    summary_file = get_nested(config, [("summary_file",), ("outputs", "summary_file")], default="outputs/qc_reports/metorigin_matching_summary.csv")
    unmatched_file = get_nested(config, [("unmatched_file",), ("outputs", "unmatched_file")], default="outputs/qc_reports/metorigin_unmatched_features.csv")
    column_report_file = get_nested(config, [("column_report_file",), ("outputs", "column_report_file")], default="outputs/qc_reports/metorigin_column_detection_report.txt")
    log_file = get_nested(config, [("log_file",), ("outputs", "log_file")], default="outputs/logs/metorigin_match.log")
    case_insensitive_names = bool(get_nested(config, [("case_insensitive_name_match",), ("matching", "case_insensitive_name_match")], default=True))

    maps_path = resolve_path(str(maps_file), project_dir)
    metorigin_path = resolve_path(str(metorigin_file), project_dir)
    output_path = resolve_path(str(output_file), project_dir)
    summary_path = resolve_path(str(summary_file), project_dir)
    unmatched_path = resolve_path(str(unmatched_file), project_dir)
    column_report_path = resolve_path(str(column_report_file), project_dir)
    log_path = resolve_path(str(log_file), project_dir)

    setup_logging(log_path)
    logging.info("Starting MAPS -> MetOriginDB matching.")
    logging.info("MAPS file: %s", maps_path)
    logging.info("MetOriginDB file: %s", metorigin_path)

    if not maps_path.exists():
        raise FileNotFoundError(f"MAPS file not found: {maps_path}")
    if not metorigin_path.exists():
        raise FileNotFoundError(f"MetOriginDB file not found: {metorigin_path}")

    maps = pd.read_csv(maps_path, low_memory=False)
    metorigin = pd.read_csv(metorigin_path, low_memory=False)
    maps = clean_dataframe_columns(maps)
    metorigin = clean_dataframe_columns(metorigin)

    logging.info("MAPS rows: %s | columns: %s", len(maps), len(maps.columns))
    logging.info("MetOriginDB rows: %s | columns: %s", len(metorigin), len(metorigin.columns))

    maps_cols = detect_maps_columns(maps)
    met_cols = detect_metorigin_columns(metorigin)
    maps_cols = apply_overrides(maps_cols, get_nested(config, [("maps_column_overrides",)], default={}) or {}, maps, "MAPS")
    met_cols = apply_overrides(met_cols, get_nested(config, [("metorigin_column_overrides",)], default={}) or {}, metorigin, "MetOriginDB")
    write_column_detection_report(column_report_path, maps, metorigin, maps_cols, met_cols)

    final_full, summary = run_matching(maps, metorigin, maps_cols, met_cols, case_insensitive_names=case_insensitive_names)
    final_full = add_origin_summary_columns(final_full)
    final_full_no_helpers = drop_internal_match_columns(final_full)
    final_compact = make_compact_output(final_full)
    final_preview = make_excel_friendly_preview(final_compact)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    unmatched_path.parent.mkdir(parents=True, exist_ok=True)

    full_path = output_path.with_name(output_path.stem + "_FULL.csv")
    compact_path = output_path.with_name(output_path.stem + "_COMPACT.csv")
    preview_path = output_path.with_name(output_path.stem + "_EXCEL_PREVIEW_truncated_details.csv")

    final_compact.to_csv(output_path, index=False)
    final_compact.to_csv(compact_path, index=False)
    final_full_no_helpers.to_csv(full_path, index=False)
    final_preview.to_csv(preview_path, index=False)
    summary.to_csv(summary_path, index=False)
    unmatched = final_compact[final_compact["metorigin_matched"] == False].copy()  # noqa: E712
    unmatched.to_csv(unmatched_path, index=False)

    logging.info("Saved main compact output: %s", output_path)
    logging.info("Saved compact output: %s", compact_path)
    logging.info("Saved full output: %s", full_path)
    logging.info("Saved Excel-friendly preview with truncated details: %s", preview_path)
    logging.info("Saved summary: %s", summary_path)
    logging.info("Saved unmatched features: %s", unmatched_path)
    logging.info("Saved column detection report: %s", column_report_path)
    logging.info("Matching summary:")
    logging.info("\n%s", summary.to_string(index=False))
    logging.info("Done.")


if __name__ == "__main__":
    main()
