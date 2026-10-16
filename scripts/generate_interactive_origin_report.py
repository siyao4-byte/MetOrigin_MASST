#!/usr/bin/env python
"""Generate an interactive HTML origin report for one or more run folders."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
from itertools import zip_longest

import pandas as pd
import yaml


ORIGIN_ORDER = ["Food", "Plant", "Human", "Tissue/Human", "Microbe", "Bacteria", "Animal", "Drug", "Environment"]
MAX_BACTERIA_LINEAGES_PER_FEATURE = 40
MAX_OTHER_ORIGIN_ITEMS_PER_FEATURE = 25


def read_csv_if_exists(path: Path, **kwargs: Any) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    return pd.read_csv(path, low_memory=False, **kwargs)


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def clean_value(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value)
    if text.lower() in {"nan", "na", "none", "unknown"}:
        return ""
    return text


def clean_smiles_key(value: Any) -> str:
    return clean_value(value).strip()


def compact_text(value: Any, limit: int = 420) -> str:
    text = clean_value(value)
    if len(text) <= limit:
        return text
    return text[: limit - 3] + "..."


def row_value(row: pd.Series, names: Iterable[str]) -> str:
    for name in names:
        if name in row.index:
            value = clean_value(row[name])
            if value:
                return value
    return ""


def load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data if isinstance(data, dict) else {}


def find_metorigin_file(run_dir: Path) -> Optional[Path]:
    candidates = [
        run_dir / "outputs" / "metorigin_matches" / "MAPS_with_MetOriginDB_origin_FULL.csv",
        run_dir / "outputs" / "metorigin_matches" / "MAPS_with_MetOriginDB_origin_COMPACT.csv",
        run_dir / "outputs" / "metorigin_matches" / "MAPS_with_MetOriginDB_origin.csv",
    ]
    found = first_existing(candidates)
    if found:
        return found
    matches = list((run_dir / "outputs").glob("**/*MetOriginDB_origin*EXCEL_PREVIEW*.csv"))
    return matches[0] if matches else None


def find_metorigin_full_file(run_dir: Path) -> Optional[Path]:
    candidates = [
        run_dir / "outputs" / "metorigin_matches" / "MAPS_with_MetOriginDB_origin_FULL.csv",
        run_dir / "outputs" / "metorigin_matches" / "MAPS_with_MetOriginDB_origin.csv",
    ]
    found = first_existing(candidates)
    if found:
        return found
    matches = list((run_dir / "outputs").glob("**/*MetOriginDB_origin_FULL*.csv"))
    return matches[0] if matches else None


def detect_analysis(run_dir: Path) -> str:
    has_metorigin = find_metorigin_file(run_dir) is not None
    has_masst = (run_dir / "Graph2_domainMASST_R_outputs").exists() or (run_dir / "Graph3_integrated_origin_R_outputs").exists()
    if has_metorigin and has_masst:
        return "MetOrigin + domainMASST"
    if has_metorigin:
        return "MetOrigin only"
    if has_masst:
        return "domainMASST only"
    return "No recognized analysis outputs"


def origin_rank(origin: str) -> int:
    try:
        return ORIGIN_ORDER.index(origin)
    except ValueError:
        return len(ORIGIN_ORDER)


def split_detail(value: Any) -> List[str]:
    text = clean_value(value)
    if not text:
        return []
    text = text.replace("{}", "|").replace("{", "|").replace("}", "|").replace(";", "|")
    text = text.replace("...", "")
    out = []
    seen = set()
    for part in text.split("|"):
        part = part.strip()
        if not part or part.lower() in {"unknown", "na", "nan", "none", "null", "not available"}:
            continue
        if part not in seen:
            seen.add(part)
            out.append(part)
    return out


def build_lineage_rows(run_dir: Path) -> List[Dict[str, str]]:
    """Read synchronized, validated lineage records; limits affect HTML only."""
    directory = run_dir / "outputs" / "metorigin_matches"
    path = directory / "bacterial_lineages.csv"
    if not path.exists():
        raise FileNotFoundError("Rebuild this run with the indexed matcher; flattened lineage columns are unsafe.")
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    fields = ["feature_id", "organism", "ncbi_id", "phylum", "class_name", "tax_order", "family", "genus", "species"]
    df = df.drop_duplicates(fields)
    rows = []
    for fid, group in df.groupby("feature_id", sort=False):
        for item in group.head(MAX_BACTERIA_LINEAGES_PER_FEATURE).to_dict("records"):
            rows.append({**item, "origin_domain": "Bacteria", "source_name": item["organism"],
                         "order": item["tax_order"], "detail": item["species"] or item["organism"]})
        if len(group) > MAX_BACTERIA_LINEAGES_PER_FEATURE:
            rows.append({"feature_id": fid, "origin_domain": "Bacteria", "source_name": "Additional complete records in bacterial_lineages.csv",
                         "detail": str(len(group) - MAX_BACTERIA_LINEAGES_PER_FEATURE), "phylum": "Additional records (display limit)"})
    associations = pd.read_csv(directory / "source_associations.csv", dtype=str, keep_default_na=False)
    associations = associations[(associations.status == "Yes") & (associations.source != "bacteria")]
    associations = associations.drop_duplicates(["feature_id", "source", "detail"])
    for (fid, source), group in associations.groupby(["feature_id", "source"], sort=False):
        for item in group.head(MAX_OTHER_ORIGIN_ITEMS_PER_FEATURE).to_dict("records"):
            rows.append({**item, "origin_domain": source.title(), "source_name": item["detail"]})
        if len(group) > MAX_OTHER_ORIGIN_ITEMS_PER_FEATURE:
            rows.append({"feature_id": fid, "origin_domain": source.title(), "source_name": "Additional records in source_associations.csv",
                         "detail": str(len(group) - MAX_OTHER_ORIGIN_ITEMS_PER_FEATURE)})
    return rows


def build_payload(run_dir: Path, project_dir: Path) -> Dict[str, Any]:
    run_name = run_dir.name
    metorigin_file = find_metorigin_file(run_dir)
    metorigin_df = read_csv_if_exists(metorigin_file) if metorigin_file else pd.DataFrame()
    if not metorigin_df.empty and "identity_id" in metorigin_df:
        flag_columns = [f"from_{s}" for s in ("human", "bacteria", "plant", "animal", "food", "drug", "environment")]
        # Aggregate before choosing a display representative; keep all source evidence.
        combined_flags = metorigin_df.groupby("identity_id")[flag_columns].agg(
            lambda s: "Yes" if s.eq("Yes").any() else "Unknown" if s.eq("Unknown").any() else "No")
        metorigin_df = metorigin_df.sort_values("metorigin_matched", ascending=False).drop_duplicates("identity_id").copy()
        metorigin_df["feature_id"] = metorigin_df["identity_id"]
        for column in flag_columns:
            metorigin_df[column] = metorigin_df.identity_id.map(combined_flags[column])
        metorigin_df["metorigin_origin_categories"] = metorigin_df.apply(
            lambda row: ";".join(column[5:] for column in flag_columns if row[column] == "Yes"), axis=1)

    graph1_long_path = run_dir / "Graph1_MetOrigin_R_outputs" / "Graph1a_feature_general_subclass_detail_long_table.csv"
    graph1_long = read_csv_if_exists(graph1_long_path)
    if graph1_long is None:
        graph1_long = pd.DataFrame()
    elif not graph1_long.empty:
        # This is an HTML display projection only. The complete CSV remains unchanged.
        graph1_long = graph1_long.drop_duplicates(["feature_id", "general_class", "layer3", "layer4"])
        graph1_long = graph1_long.groupby(["feature_id", "general_class"], sort=False).head(MAX_OTHER_ORIGIN_ITEMS_PER_FEATURE)

    graph2_long_path = run_dir / "Graph2_domainMASST_R_outputs" / "Graph2b_domainMASST_sankey_long_table.csv"
    graph2_long = read_csv_if_exists(graph2_long_path)
    if graph2_long is None:
        graph2_long = pd.DataFrame()

    graph3_evidence_path = run_dir / "Graph3_integrated_origin_R_outputs" / "Graph3_integrated_origin_evidence_long.csv"
    graph3_evidence = read_csv_if_exists(graph3_evidence_path)
    if graph3_evidence is None:
        graph3_evidence = pd.DataFrame()
    graph3_sankey_path = run_dir / "Graph3_integrated_origin_R_outputs" / "Graph3c_integrated_sankey_long_table.csv"
    graph3_sankey = read_csv_if_exists(graph3_sankey_path)
    if graph3_sankey is None:
        graph3_sankey = pd.DataFrame()
    lineage_rows = build_lineage_rows(run_dir)

    summary_path = run_dir / "outputs" / "qc_reports" / "metorigin_matching_summary.csv"
    summary_df = read_csv_if_exists(summary_path)
    if summary_df is None:
        summary_df = pd.DataFrame()

    root_config = load_yaml(project_dir / "config" / "config.yaml")
    run_config_files = list((run_dir / "config").glob("*.yaml"))
    run_config = load_yaml(run_config_files[0]) if run_config_files else {}
    config = {**root_config, **run_config}
    np_cache_path = project_dir / "data" / "npclassifier_cache.csv"
    np_cache: Dict[str, Dict[str, str]] = {}
    if np_cache_path.exists():
        np_df = pd.read_csv(np_cache_path, low_memory=False)
        if "SMILES" in np_df.columns:
            for _, cache_row in np_df.dropna(subset=["SMILES"]).iterrows():
                key = clean_smiles_key(cache_row.get("SMILES"))
                if not key or key in np_cache:
                    continue
                np_cache[key] = {
                    "np_pathway": clean_value(cache_row.get("np_pathway")),
                    "np_superclass": clean_value(cache_row.get("np_superclass")),
                    "np_class": clean_value(cache_row.get("np_class")),
                    "is_glycoside": clean_value(cache_row.get("is_glycoside")),
                }

    features = []
    if not metorigin_df.empty:
        metorigin_df["feature_id"] = metorigin_df["feature_id"].astype(str)
        for _, row in metorigin_df.iterrows():
            feature_id = clean_value(row.get("feature_id"))
            compound_class = row_value(
                row,
                [
                    "canopus_npc_pathway",
                    "canopus_npc_superclass",
                    "canopus_classyfire_superclass",
                    "canopus_classyfire_subclass",
                ],
            )
            smiles = row_value(row, ["smiles", "metorigin_smiles_id"])
            np_hit = {key: row_value(row, [key]) for key in ("np_pathway", "np_superclass", "np_class")}
            cached = np_cache.get(clean_smiles_key(smiles), {})
            np_hit = {key: value or cached.get(key, "") for key, value in np_hit.items()}
            np_compound_class = row_value(pd.Series(np_hit), ["np_pathway", "np_superclass", "np_class"]) if np_hit else ""
            compound_class = np_compound_class or compound_class
            feature_label = row_value(row, ["compound_name", "feature_label", "feature_usi"])
            if not feature_label:
                mz = clean_value(row.get("mz"))
                rt = clean_value(row.get("rt"))
                feature_label = f"{mz}_{rt}_{feature_id}" if mz or rt else feature_id
            features.append(
                {
                    "feature_id": feature_id,
                    "feature_label": f"{feature_label} | {feature_id}",
                    "compound_name": row_value(row, ["compound_name", "metorigin_compound_name"]),
                    "smiles": compact_text(smiles, 650),
                    "cid": row_value(row, ["cid", "metorigin_pubchem_compound_id"]),
                    "hmdb_id": row_value(row, ["hmdb_id", "metorigin_hmdb_id"]),
                    "kegg_id": row_value(row, ["metorigin_kegg_id", "kegg_id"]),
                    "annotation_type": row_value(row, ["annotation_type", "propagated_annotation_type"]),
                    "annotation_source": row_value(row, ["metorigin_match_method", "annotation_type"]) or "MetOrigin",
                    "confidence_level": row_value(row, ["confidence_level"]),
                    "confidence_score": row_value(row, ["confidence_score"]),
                    "compound_class": compound_class,
                    "np_pathway": np_hit.get("np_pathway", ""),
                    "np_superclass": np_hit.get("np_superclass", ""),
                    "np_class": np_hit.get("np_class", ""),
                    "np_is_glycoside": np_hit.get("is_glycoside", ""),
                    "superclass": np_hit.get("np_superclass", "") or row_value(row, ["canopus_npc_superclass", "canopus_classyfire_superclass"]),
                    "subclass": np_hit.get("np_class", "") or row_value(row, ["canopus_classyfire_subclass", "canopus_classyfire_specclass"]),
                    "metorigin_matched": row_value(row, ["metorigin_matched"]),
                    "metorigin_origin_categories": row_value(row, ["metorigin_origin_categories"]),
                    "metorigin_n_database_records": row_value(row, ["metorigin_n_database_records"]),
                    "from_food": row_value(row, ["from_food"]),
                    "from_human": row_value(row, ["from_human"]),
                    "from_bacteria": row_value(row, ["from_bacteria"]),
                    "from_plant": row_value(row, ["from_plant"]),
                    "from_animal": row_value(row, ["from_animal"]),
                    "from_drug": row_value(row, ["from_drug"]),
                    "from_environment": row_value(row, ["from_environment"]),
                    "food_detail": compact_text(row_value(row, ["from_which_food"]), 800),
                    "human_detail": compact_text(row_value(row, ["from_which_part"]), 800),
                    "microbiome_detail": compact_text(row_value(row, ["bacteria_genus", "from_which_bacteria", "bacteria_species"]), 800),
                    "plant_detail": compact_text(row_value(row, ["from_which_plant"]), 800),
                    "animal_detail": compact_text(row_value(row, ["from_which_animal"]), 800),
                    "drug_detail": compact_text(row_value(row, ["from_which_drug"]), 800),
                    "environment_detail": compact_text(row_value(row, ["from_which_environment"]), 800),
                }
            )

    feature_ids_from_metorigin = {f["feature_id"] for f in features}

    origins = []
    if not graph1_long.empty:
        graph1_long["feature_id"] = graph1_long["feature_id"].astype(str)
        for _, row in graph1_long.iterrows():
            feature_id = clean_value(row.get("feature_id"))
            origins.append(
                {
                    "feature_id": feature_id,
                    "feature_label": row_value(row, ["feature_label"]),
                    "source": "MetOrigin",
                    "general_class": row_value(row, ["general_class"]),
                    "subclass": row_value(row, ["layer3_plot", "layer3"]),
                    "detail": row_value(row, ["layer4_plot", "layer4"]),
                    "confidence_level": row_value(row, ["confidence_level"]),
                    "annotation_source": row_value(row, ["annotation_type"]),
                }
            )

    masst_rows = []
    if not graph2_long.empty:
        graph2_long["feature_id"] = graph2_long["feature_id"].astype(str)
        for _, row in graph2_long.iterrows():
            feature_id = clean_value(row.get("feature_id"))
            domain = row_value(row, ["domain_label", "domain"]).title()
            masst_rows.append(
                {
                    "feature_id": feature_id,
                    "source": "domainMASST",
                    "domain": domain,
                    "subclass": row_value(row, ["subclass_name"]),
                    "source_label": row_value(row, ["source_label_plot", "source_label"]),
                    "evidence_weight": row_value(row, ["evidence_weight"]),
                    "occurrence_fraction": row_value(row, ["max_occurrence_fraction"]),
                }
            )
            if feature_id not in feature_ids_from_metorigin:
                label = row_value(row, ["feature_label", "feature_label_raw", "mz_rt_feature_label"]) or feature_id
                features.append(
                    {
                        "feature_id": feature_id,
                        "feature_label": label,
                        "compound_name": label.split(" | ")[0],
                        "smiles": "",
                        "cid": "",
                        "hmdb_id": "",
                        "kegg_id": "",
                        "annotation_type": "domainMASST",
                        "annotation_source": "domainMASST",
                        "confidence_level": "",
                        "confidence_score": "",
                        "compound_class": "",
                        "np_pathway": "",
                        "np_superclass": "",
                        "np_class": "",
                        "np_is_glycoside": "",
                        "superclass": "",
                        "subclass": "",
                        "metorigin_matched": "FALSE",
                        "metorigin_origin_categories": "",
                        "metorigin_n_database_records": "",
                    }
                )
                feature_ids_from_metorigin.add(feature_id)

    evidence_rows = []
    if not graph3_sankey.empty:
        graph3_sankey["feature_id"] = graph3_sankey["feature_id"].astype(str)
        for _, row in graph3_sankey.iterrows():
            evidence_rows.append(
                {
                    "feature_id": clean_value(row.get("feature_id")),
                    "feature_label": row_value(row, ["feature_label"]),
                    "evidence_source": row_value(row, ["evidence_source"]),
                    "origin_domain": row_value(row, ["origin_domain_normalized"]),
                    "subclass": row_value(row, ["source_label_for_plot"]),
                    "detail": row_value(row, ["detail_label_plot", "source_label_for_plot"]),
                    "weight": "1",
                }
            )
    else:
        for row in origins:
            evidence_rows.append(
                {
                    "feature_id": row["feature_id"],
                    "feature_label": row.get("feature_label", ""),
                    "evidence_source": "MetOrigin",
                    "origin_domain": row.get("general_class", ""),
                    "subclass": row.get("subclass", ""),
                    "detail": row.get("detail", ""),
                    "weight": "1",
                }
            )
        for row in masst_rows:
            evidence_rows.append(
                {
                    "feature_id": row["feature_id"],
                    "feature_label": "",
                    "evidence_source": "domainMASST",
                    "origin_domain": row.get("domain", ""),
                    "subclass": row.get("subclass", ""),
                    "detail": row.get("source_label", ""),
                    "weight": row.get("evidence_weight", "1") or "1",
                }
            )

    masst_summary: Dict[str, Dict[str, Any]] = {}
    for row in masst_rows:
        fid = row["feature_id"]
        entry = masst_summary.setdefault(fid, {"domains": set(), "subclasses": set(), "microbe_hits": 0, "human_plasma_hits": 0})
        entry["domains"].add(row["domain"])
        entry["subclasses"].add(row["subclass"])
        try:
            weight = int(float(row["evidence_weight"]))
        except Exception:
            weight = 0
        if row["domain"].lower() in {"microbe", "bacteria"}:
            entry["microbe_hits"] += weight
        if "plasma" in row["source_label"].lower() or row["domain"].lower() in {"human", "tissue/human"}:
            entry["human_plasma_hits"] += weight

    for entry in masst_summary.values():
        entry["domains"] = "; ".join(sorted(entry["domains"]))
        entry["subclasses"] = "; ".join(sorted(s for s in entry["subclasses"] if s))

    if graph3_evidence.empty:
        integrated_count = 0
    else:
        integrated_count = int(graph3_evidence["feature_id"].astype(str).nunique())

    feature_count = len({f["feature_id"] for f in features})
    metorigin_matched_count = len({f["feature_id"] for f in features if str(f.get("metorigin_matched", "")).lower() == "true"})
    masst_feature_count = len(masst_summary)
    no_hit_count = max(feature_count - len({*masst_summary.keys(), *[f["feature_id"] for f in features if str(f.get("metorigin_matched", "")).lower() == "true"]}), 0)

    domain_params = config.get("domainmasst", {}) if isinstance(config.get("domainmasst", {}), dict) else {}
    metorigin_params = {
        "maps_annotation_file": config.get("maps_annotation_file", ""),
        "metorigin_file": config.get("metorigin_file", ""),
        "matching_priority": config.get("matching_priority", ""),
        "case_insensitive_name_match": config.get("case_insensitive_name_match", config.get("case_insensitive_name_matching", "")),
        "exact_smiles_only": config.get("exact_smiles_only", ""),
    }

    return {
        "meta": {
            "run_name": run_name,
            "run_dir": str(run_dir),
            "dataset": run_name,
            "analysis": detect_analysis(run_dir),
            "metorigin_file": str(metorigin_file or ""),
            "graph1_long_table": str(graph1_long_path if graph1_long_path.exists() else ""),
            "domainmasst_table": str(graph2_long_path if graph2_long_path.exists() else ""),
            "integrated_evidence_table": str(graph3_evidence_path if graph3_evidence_path.exists() else ""),
            "integrated_sankey_table": str(graph3_sankey_path if graph3_sankey_path.exists() else ""),
        },
        "params": {
            "metorigin": {**metorigin_params, "npclassifier_cache": str(np_cache_path if np_cache_path.exists() else "")},
            "domainmasst": domain_params,
        },
        "summary": {
            "feature_count": feature_count,
            "metorigin_matched_features": metorigin_matched_count,
            "domainmasst_matched_features": masst_feature_count,
            "integrated_evidence_features": integrated_count,
            "features_without_any_hit": no_hit_count,
            "metorigin_summary_rows": summary_df.fillna("").to_dict(orient="records") if not summary_df.empty else [],
        },
        "features": features,
        "origins": origins,
        "evidence": evidence_rows,
        "masst": masst_rows,
        "masst_summary": masst_summary,
        "lineage": lineage_rows,
    }


def write_report(payload: Dict[str, Any], out_path: Path) -> None:
    data_json = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    title = html.escape(f"Recorded Source Associations - {payload['meta']['run_name']}")
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root { --ink:#1f2933; --muted:#667085; --line:#d0d5dd; --bg:#f6f8fb; --panel:#ffffff; }
    * { box-sizing: border-box; }
    body { margin:0; font-family: Arial, Helvetica, sans-serif; color:var(--ink); background:var(--bg); }
    header { padding:24px 32px 18px; background:#fff; border-bottom:1px solid var(--line); }
    h1 { margin:0 0 8px; font-size:26px; }
    h2 { margin:0 0 14px; font-size:18px; }
    .meta-grid, .controls { display:grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap:12px; }
    .metric { border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:12px; }
    .metric b { display:block; font-size:12px; color:var(--muted); text-transform:uppercase; margin-bottom:4px; }
    main { padding:20px 32px 32px; }
    section { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:16px; margin-bottom:18px; }
    label { display:block; font-size:12px; font-weight:700; color:#475467; margin-bottom:4px; }
    select, input { width:100%; border:1px solid var(--line); border-radius:6px; padding:8px 10px; background:#fff; color:var(--ink); }
    .checkbox-group { min-height:38px; max-height:116px; overflow:auto; border:1px solid var(--line); border-radius:6px; padding:6px 8px; background:#fff; }
    .checkbox-group label { display:inline-flex; align-items:center; gap:5px; margin:3px 10px 3px 0; font-size:12px; font-weight:400; color:var(--ink); }
    .checkbox-group input { width:auto; }
    .plot { width:100%; height:720px; }
    #sankeyPlot { height:1080px; }
    .split { display:grid; grid-template-columns: minmax(0, 1fr) 720px; gap:16px; align-items:start; }
    #detailPanel { position:sticky; top:12px; max-height:720px; overflow:auto; }
    .hint { color:var(--muted); font-size:13px; }
    table { width:100%; border-collapse:collapse; font-size:13px; }
    th, td { border-bottom:1px solid #eaecf0; padding:7px 8px; text-align:left; vertical-align:top; }
    th { background:#f9fafb; font-weight:700; }
    .pill { display:inline-block; padding:2px 7px; border-radius:999px; background:#eef2ff; margin:2px 3px 2px 0; font-size:12px; }
    .small { font-size:12px; color:var(--muted); }
    .lineage-tree { max-height:320px; overflow:auto; border:1px solid #eaecf0; border-radius:6px; padding:8px 10px; background:#fcfcfd; font-size:12px; }
    .lineage-tree ul { margin:4px 0 4px 16px; padding:0; }
    .lineage-tree li { margin:2px 0; }
    @media (max-width: 1100px) { .split { grid-template-columns: 1fr; } #detailPanel { position:static; } }
  </style>
</head>
<body>
<header>
  <h1>Interactive Origin Report</h1>
  <div class="meta-grid">
    <div class="metric"><b>Dataset</b><span id="datasetName"></span></div>
    <div class="metric"><b>Analysis Run</b><span id="analysisName"></span></div>
    <div class="metric"><b>Features</b><span id="featureCount"></span></div>
    <div class="metric"><b>Without Any Hit</b><span id="noHitCount"></span></div>
  </div>
</header>
<main>
  <section>
    <h2>Run Parameters</h2>
    <div class="split">
      <div><h3>MetOrigin</h3><table id="metParams"></table></div>
      <div><h3>domainMASST</h3><table id="masstParams"></table></div>
    </div>
  </section>
  <section>
    <h2>Filters</h2>
    <div class="controls">
      <div><label for="confidenceFilter">Confidence level</label><select id="confidenceFilter"></select></div>
      <div><label for="evidenceFilter">Origin evidence</label><select id="evidenceFilter"></select></div>
      <div><label>Recorded source</label><div id="originFilter" class="checkbox-group"></div></div>
      <div><label for="sourceFilter">Annotation source</label><select id="sourceFilter"></select></div>
      <div><label for="classFilter">NPClassifier class</label><select id="classFilter"></select></div>
      <div><label for="searchBox">Search compound / SMILES / CID / ID</label><input id="searchBox" placeholder="type any text"></div>
    </div>
    <p class="hint" id="filterNote"></p>
  </section>
  <section>
    <h2>Interactive Sankey</h2>
    <p class="hint">Recorded source associations, not definitive biological origins. Human-only means only human evidence is recorded under the matching criteria, not exclusive human production. Click a feature for details.</p>
    <p class="hint">Display limits never change analytical counts. Complete records: <a href="outputs/metorigin_matches/source_associations.csv">source associations</a> and <a href="outputs/metorigin_matches/bacterial_lineages.csv">bacterial lineages</a>.</p>
    <div id="sankeyPlot" class="plot"></div>
  </section>
  <section>
    <h2>Interactive Circular Network</h2>
    <p class="hint">Outer features are sorted by primary origin. Center-node size is proportional to unique feature count after filtering.</p>
    <div class="split">
      <div id="circularPlot" class="plot"></div>
      <aside id="detailPanel" class="metric"><b>Clicked feature</b><div id="featureDetails" class="hint">Click a feature node in either plot.</div></aside>
    </div>
  </section>
  <section>
    <h2>Summary Table</h2>
    <table id="summaryTable"></table>
  </section>
</main>
<script>
const DATA = __DATA__;
const ORIGIN_COLORS = {Food:"#E9C75F", Plant:"#78BE7C", Human:"#DC7C78", "Tissue/Human":"#DC7C78", Microbe:"#5DA7D7", Bacteria:"#5DA7D7", Animal:"#B99A78", Drug:"#A678D6", Environment:"#87969B", domainMASST:"#8EC5EA", MetOrigin:"#1F5F91"};
const ORIGIN_ORDER = ["Food","Plant","Human","Tissue/Human","Microbe","Bacteria","Animal","Drug","Environment"];
const MAX_SANKEY_FEATURES = 350;

const featureMap = new Map(DATA.features.map(f => [String(f.feature_id), f]));
const originByFeature = new Map();
for (const o of DATA.origins) {
  const fid = String(o.feature_id);
  if (!originByFeature.has(fid)) originByFeature.set(fid, []);
  originByFeature.get(fid).push(o);
}
const evidenceByFeature = new Map();
for (const e of DATA.evidence || []) {
  const fid = String(e.feature_id);
  if (!evidenceByFeature.has(fid)) evidenceByFeature.set(fid, []);
  evidenceByFeature.get(fid).push(e);
}
const masstByFeature = new Map();
for (const m of DATA.masst) {
  const fid = String(m.feature_id);
  if (!masstByFeature.has(fid)) masstByFeature.set(fid, []);
  masstByFeature.get(fid).push(m);
}
const lineageByFeature = new Map();
for (const l of DATA.lineage || []) {
  const fid = String(l.feature_id);
  if (!lineageByFeature.has(fid)) lineageByFeature.set(fid, []);
  lineageByFeature.get(fid).push(l);
}

function esc(x) {
  return String(x ?? "").replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
}
function uniq(values) {
  return [...new Set(values.filter(v => v !== undefined && v !== null && String(v).trim() !== ""))].sort((a,b) => String(a).localeCompare(String(b), undefined, {numeric:true}));
}
function tableFromObject(obj) {
  const rows = Object.entries(obj || {}).map(([k,v]) => `<tr><th>${esc(k)}</th><td>${esc(Array.isArray(v) ? v.join("; ") : JSON.stringify(v) === "{}" ? "" : v)}</td></tr>`).join("");
  return rows || "<tr><td class='hint'>No parameters detected</td></tr>";
}
function fillSelect(id, values, allLabel) {
  const el = document.getElementById(id);
  el.innerHTML = `<option value="">${allLabel}</option>` + values.map(v => `<option value="${esc(v)}">${esc(v)}</option>`).join("");
}
function fillCheckboxGroup(id, values) {
  const el = document.getElementById(id);
  el.innerHTML = values.map(v => `<label><input type="checkbox" value="${esc(v)}"> ${esc(v)}</label>`).join("");
}
function primaryOrigin(fid) {
  const evidence = selectedEvidenceForFeature(fid);
  if (!evidence.length) return "";
  return evidence.map(o => o.origin_domain).sort((a,b) => ORIGIN_ORDER.indexOf(a) - ORIGIN_ORDER.indexOf(b))[0] || evidence[0].origin_domain;
}
function selectedEvidenceSource() {
  const el = document.getElementById("evidenceFilter");
  return el ? el.value : "";
}
function selectedOrigins() {
  return [...document.querySelectorAll("#originFilter input:checked")].map(x => x.value);
}
function selectedEvidenceForFeature(fid) {
  const source = selectedEvidenceSource();
  const origins = selectedOrigins();
  let rows = evidenceByFeature.get(String(fid)) || [];
  if (source) rows = rows.filter(r => r.evidence_source === source);
  if (origins.length) rows = rows.filter(r => origins.includes(r.origin_domain));
  return rows;
}
function filteredFeatures() {
  const conf = document.getElementById("confidenceFilter").value;
  const src = document.getElementById("sourceFilter").value;
  const cls = document.getElementById("classFilter").value;
  const q = document.getElementById("searchBox").value.trim().toLowerCase();
  return DATA.features.filter(f => {
    if (conf && String(f.confidence_level) !== conf) return false;
    if (src && String(f.annotation_source) !== src && String(f.annotation_type) !== src) return false;
    if (cls && ![f.np_pathway, f.np_superclass, f.np_class].map(String).includes(cls)) return false;
    if (q) {
      const hay = [f.feature_id, f.feature_label, f.compound_name, f.smiles, f.cid, f.hmdb_id, f.kegg_id].join(" ").toLowerCase();
      if (!hay.includes(q)) return false;
    }
    if ((selectedEvidenceSource() || selectedOrigins().length) && !selectedEvidenceForFeature(f.feature_id).length) return false;
    if (!selectedEvidenceSource() && !selectedOrigins().length && DATA.evidence && DATA.evidence.length && !selectedEvidenceForFeature(f.feature_id).length) return false;
    return true;
  }).sort((a,b) => {
    const ao = primaryOrigin(a.feature_id), bo = primaryOrigin(b.feature_id);
    const ar = ORIGIN_ORDER.indexOf(ao); const br = ORIGIN_ORDER.indexOf(bo);
    return (ar < 0 ? 99 : ar) - (br < 0 ? 99 : br) || String(a.feature_label).localeCompare(String(b.feature_label));
  });
}
function featureHover(f) {
  return [
    `<b>${esc(f.compound_name || f.feature_label)}</b>`,
    `Feature ID: ${esc(f.feature_id)}`,
    `Confidence: ${esc(f.confidence_level)}`,
    `Annotation source: ${esc(f.annotation_source)}`,
    `NP pathway: ${esc(f.np_pathway)}`,
    `NP superclass: ${esc(f.np_superclass)}`,
    `NP class: ${esc(f.np_class)}`,
    `Superclass: ${esc(f.superclass)}`,
    `Subclass: ${esc(f.subclass)}`,
    `SMILES: ${esc(f.smiles)}`,
    `CID: ${esc(f.cid)}`,
    `HMDB: ${esc(f.hmdb_id)}`
  ].join("<br>");
}
function renderSankey(features) {
  const selected = features.slice(0, MAX_SANKEY_FEATURES);
  const nodes = [], labels = [], colors = [], custom = [];
  const nodeIndex = new Map();
  function addNode(key, label, color, meta) {
    if (nodeIndex.has(key)) return nodeIndex.get(key);
    const idx = nodes.length; nodeIndex.set(key, idx); nodes.push(key); labels.push(label); colors.push(color); custom.push(meta); return idx;
  }
  function addNodeFeature(idx, fid) {
    if (!custom[idx].feature_ids) custom[idx].feature_ids = [];
    if (!custom[idx].feature_ids.includes(String(fid))) custom[idx].feature_ids.push(String(fid));
  }
  const links = new Map();
  function addLink(s,t,val,color) {
    const key = `${s}|${t}|${color}`;
    if (!links.has(key)) links.set(key, {source:s, target:t, value:0, color});
    links.get(key).value += val;
  }
  for (const f of selected) {
    const sourceRows = selectedEvidenceForFeature(f.feature_id);
    const fidx = addNode(`F:${f.feature_id}`, f.compound_name || f.feature_label, "#3A3A3A", {type:"feature", feature_id:String(f.feature_id), feature_ids:[String(f.feature_id)], hover:featureHover(f)});
    for (const o of sourceRows) {
      const domain = o.origin_domain || "Origin not specified";
      const color = ORIGIN_COLORS[domain] || (o.evidence_source === "domainMASST" ? ORIGIN_COLORS.domainMASST : "#888");
      const sourceDomain = selectedEvidenceSource() ? domain : `${o.evidence_source}: ${domain}`;
      const oidx = addNode(`O:${o.evidence_source}:${domain}`, sourceDomain, color, {type:"origin", evidence_source:o.evidence_source, origin_domain:domain, feature_ids:[], hover:`<b>${esc(o.evidence_source)}</b><br>${esc(domain)}`});
      const detailLabel = o.detail || o.subclass || "Unspecified detail";
      const sidx = addNode(`S:${o.evidence_source}:${domain}:${detailLabel}`, detailLabel, color, {type:"subclass", evidence_source:o.evidence_source, origin_domain:domain, detail:detailLabel, feature_ids:[], hover:`${esc(o.evidence_source)} / ${esc(domain)} / ${esc(detailLabel)}`});
      addNodeFeature(oidx, f.feature_id);
      addNodeFeature(sidx, f.feature_id);
      addLink(fidx, oidx, 1, color.replace(")", ",0.25)").replace("rgb", "rgba"));
      addLink(oidx, sidx, 1, color.replace(")", ",0.35)").replace("rgb", "rgba"));
    }
  }
  const linkArr = [...links.values()];
  Plotly.react("sankeyPlot", [{
    type: "sankey",
    arrangement: "snap",
    node: {label: labels, color: colors, pad: 10, thickness: 12, customdata: custom, hovertemplate: "%{customdata.hover}<extra></extra>"},
    link: {source: linkArr.map(l => l.source), target: linkArr.map(l => l.target), value: linkArr.map(l => l.value), color: linkArr.map(l => l.color), hovertemplate: "Features: %{value}<extra></extra>"}
  }], {margin:{l:10,r:10,t:10,b:10}, font:{size:10}}, {responsive:true});
  const note = selected.length < features.length ? `Showing ${selected.length} of ${features.length} filtered features in Sankey. Narrow filters or search to inspect more.` : `Showing ${features.length} filtered features.`;
  document.getElementById("filterNote").textContent = note;
  document.getElementById("sankeyPlot").on("plotly_click", ev => {
    const p = ev.points && ev.points[0];
    if (!p) return;
    const meta = p.customdata || custom[p.pointNumber];
    if (!meta) return;
    if (meta.type === "feature") showFeature(meta.feature_id);
    else showNode(meta);
  });
}
function renderCircular(features) {
  const originCounts = {};
  const ids = new Set(features.map(f => String(f.feature_id)));
  for (const f of features) {
    for (const o of selectedEvidenceForFeature(f.feature_id)) {
      const key = selectedEvidenceSource() ? o.origin_domain : `${o.evidence_source}: ${o.origin_domain}`;
      originCounts[key] = (originCounts[key] || 0) + 1;
    }
  }
  const n = features.length || 1;
  const featureX=[], featureY=[], featureText=[], featureColor=[], featureCustom=[];
  features.forEach((f,i) => {
    const a = 2 * Math.PI * i / n - Math.PI / 2;
    featureX.push(Math.cos(a)); featureY.push(Math.sin(a)); featureText.push(featureHover(f)); featureColor.push(ORIGIN_COLORS[primaryOrigin(f.feature_id)] || "#333"); featureCustom.push(String(f.feature_id));
  });
  const traces = [];
  const originKeys = Object.keys(originCounts);
  for (const origin of originKeys) {
    const domain = origin.includes(": ") ? origin.split(": ").slice(1).join(": ") : origin;
    const angle = (2 * Math.PI * (originKeys.indexOf(origin)) / Math.max(originKeys.length, 1)) - Math.PI / 2;
    const ox = 0.42 * Math.cos(angle), oy = 0.42 * Math.sin(angle);
    const xs=[], ys=[];
    features.forEach((f,i) => {
      if (selectedEvidenceForFeature(f.feature_id).some(o => {
        const key = selectedEvidenceSource() ? o.origin_domain : `${o.evidence_source}: ${o.origin_domain}`;
        return key === origin;
      })) {
        xs.push(featureX[i], ox, null); ys.push(featureY[i], oy, null);
      }
    });
    traces.push({type:"scattergl", mode:"lines", x:xs, y:ys, line:{color:ORIGIN_COLORS[domain] || "#888", width:0.55}, opacity:0.22, hoverinfo:"skip", showlegend:false});
  }
  traces.push({type:"scattergl", mode:"markers", x:featureX, y:featureY, marker:{size:4, color:featureColor, opacity:0.9}, text:featureText, customdata:featureCustom, hovertemplate:"%{text}<extra></extra>", name:"Features", showlegend:false});
  const ox=[], oy=[], osize=[], otext=[], ocolor=[];
  for (const origin of originKeys) {
    const domain = origin.includes(": ") ? origin.split(": ").slice(1).join(": ") : origin;
    const angle = (2 * Math.PI * (originKeys.indexOf(origin)) / Math.max(originKeys.length, 1)) - Math.PI / 2;
    ox.push(0.42 * Math.cos(angle)); oy.push(0.42 * Math.sin(angle)); osize.push(18 + Math.sqrt(originCounts[origin]) * 0.8); otext.push(`<b>${origin}</b><br>Feature count: ${originCounts[origin]}`); ocolor.push(ORIGIN_COLORS[domain] || "#888");
  }
  traces.push({
    type:"scatter",
    mode:"markers",
    x:ox,
    y:oy,
    marker:{size:osize.map(s => s + 12), color:"rgba(255,255,255,0.92)", line:{color:"rgba(30,30,30,0.28)", width:1.5}},
    hoverinfo:"skip",
    name:"Origin node halo",
    showlegend:false
  });
  traces.push({
    type:"scatter",
    mode:"markers+text",
    x:ox,
    y:oy,
    marker:{size:osize.map(s => s + 4), color:ocolor, opacity:0.96, line:{color:"#111111", width:2.2}},
    text:originKeys,
    textfont:{size:14, color:"#111111"},
    textposition:"middle right",
    hovertext:otext,
    hovertemplate:"%{hovertext}<extra></extra>",
    name:"Origin nodes",
    showlegend:false
  });
  Plotly.react("circularPlot", traces, {margin:{l:10,r:10,t:10,b:10}, xaxis:{visible:false, range:[-1.18,1.18]}, yaxis:{visible:false, range:[-1.18,1.18], scaleanchor:"x"}, plot_bgcolor:"#fff", paper_bgcolor:"#fff"}, {responsive:true});
  document.getElementById("circularPlot").on("plotly_click", ev => {
    const p = ev.points && ev.points[0];
    if (p && p.customdata) showFeature(String(p.customdata));
  });
}
function showFeature(fid) {
  const f = featureMap.get(String(fid));
  if (!f) return;
  const origins = originByFeature.get(String(fid)) || [];
  const masst = masstByFeature.get(String(fid)) || [];
  const lineage = lineageByFeature.get(String(fid)) || [];
  const ms = DATA.masst_summary[String(fid)] || {};
  const originRows = origins.map(o => `<tr><td>${esc(o.general_class)}</td><td>${esc(o.subclass)}</td><td>${esc(o.detail)}</td></tr>`).join("") || "<tr><td colspan='3'>No MetOrigin origin rows</td></tr>";
  const masstRows = masst.map(m => `<tr><td>${esc(m.domain)}</td><td>${esc(m.subclass)}</td><td>${esc(m.evidence_weight)}</td><td>${esc(m.occurrence_fraction)}</td></tr>`).join("") || "<tr><td colspan='4'>No domainMASST rows</td></tr>";
  const lineageTree = renderLineageTree(lineage);
  document.getElementById("featureDetails").innerHTML = `
    <h3>${esc(f.compound_name || f.feature_label)}</h3>
    <p class="small">Feature ID: ${esc(f.feature_id)}</p>
    <p>${["CID: "+esc(f.cid), "HMDB: "+esc(f.hmdb_id), "KEGG: "+esc(f.kegg_id), "Confidence: "+esc(f.confidence_level), "Source: "+esc(f.annotation_source)].map(x => `<span class="pill">${x}</span>`).join("")}</p>
    <table>${tableFromObject({"NP pathway":f.np_pathway, "NP superclass":f.np_superclass, "NP class":f.np_class, "NP glycoside":f.np_is_glycoside, SMILES:f.smiles, Superclass:f.superclass, Subclass:f.subclass, "MetOrigin records":f.metorigin_n_database_records, Food:f.food_detail, Microbiome:f.microbiome_detail, "Human/plasma":f.human_detail, Plant:f.plant_detail, Animal:f.animal_detail, Drug:f.drug_detail, Environment:f.environment_detail, "domainMASST domains":ms.domains || "", "domainMASST subclasses":ms.subclasses || "", "Microbiome hits":ms.microbe_hits || 0, "Human/plasma hits":ms.human_plasma_hits || 0})}</table>
    <h4>Recorded bacterial taxonomy</h4>${lineageTree}
    <h4>Recorded source associations</h4><table><tr><th>Source</th><th>Subclass</th><th>Detail</th></tr>${originRows}</table>
    <h4>domainMASST matches</h4><table><tr><th>Domain</th><th>Phylum/Subclass</th><th>Hits</th><th>Occurrence fraction</th></tr>${masstRows}</table>`;
}
function showNode(meta) {
  const featureIds = [...new Set((meta.feature_ids || []).map(String))];
  const titleBits = [meta.evidence_source, meta.origin_domain, meta.detail].filter(Boolean);
  const scopeTitle = titleBits.join(" / ") || "Selected Sankey node";
  const features = featureIds.map(fid => featureMap.get(fid)).filter(Boolean);
  const scopedEvidence = [];
  for (const fid of featureIds) {
    for (const e of selectedEvidenceForFeature(fid)) {
      if (meta.evidence_source && e.evidence_source !== meta.evidence_source) continue;
      if (meta.origin_domain && e.origin_domain !== meta.origin_domain) continue;
      if (meta.detail && (e.detail || e.subclass || "Unspecified detail") !== meta.detail) continue;
      scopedEvidence.push(e);
    }
  }
  const subgroupCounts = {};
  for (const e of scopedEvidence) {
    const key = e.detail || e.subclass || e.origin_domain || "Unspecified";
    if (!subgroupCounts[key]) subgroupCounts[key] = new Set();
    subgroupCounts[key].add(String(e.feature_id));
  }
  const subgroupRows = Object.entries(subgroupCounts)
    .sort((a,b) => b[1].size - a[1].size || a[0].localeCompare(b[0]))
    .slice(0, 80)
    .map(([k,v]) => `<tr><td>${esc(k)}</td><td>${v.size}</td></tr>`)
    .join("") || "<tr><td colspan='2'>No subgroup rows available</td></tr>";
  const featureRows = features
    .slice(0, 160)
    .map(f => `<tr><td><button type="button" onclick="showFeature('${esc(f.feature_id)}')">open</button></td><td>${esc(f.feature_id)}</td><td>${esc(f.compound_name || f.feature_label)}</td><td>${esc(f.confidence_level)}</td><td>${esc(f.annotation_source)}</td><td>${esc(f.superclass || f.compound_class)}</td></tr>`)
    .join("") || "<tr><td colspan='6'>No connected features available</td></tr>";
  const omitted = features.length > 160 ? `<p class="hint">${features.length - 160} additional features omitted from this panel. Use filters/search to narrow this node.</p>` : "";
  document.getElementById("featureDetails").innerHTML = `
    <h3>${esc(scopeTitle)}</h3>
    <p><span class="pill">${featureIds.length} connected features</span><span class="pill">${scopedEvidence.length} evidence rows</span></p>
    <h4>Subgroups</h4>
    <table><tr><th>Subgroup/detail</th><th>Feature count</th></tr>${subgroupRows}</table>
    <h4>Connected features</h4>
    ${omitted}
    <table><tr><th></th><th>Feature ID</th><th>Compound</th><th>Confidence</th><th>Annotation source</th><th>Compound class</th></tr>${featureRows}</table>`;
}
function nestedInsert(root, path) {
  let node = root;
  for (const part of path.filter(Boolean)) {
    if (!node[part]) node[part] = {};
    node = node[part];
  }
}
function treeToHtml(obj) {
  const keys = Object.keys(obj).sort((a,b) => a.localeCompare(b, undefined, {numeric:true}));
  if (!keys.length) return "";
  return `<ul>${keys.map(k => `<li>${esc(k)}${treeToHtml(obj[k])}</li>`).join("")}</ul>`;
}
function renderLineageTree(rows) {
  if (!rows.length) return "<p class='hint'>No expanded lineage rows available for this feature.</p>";
  const bacteriaRows = rows.filter(r => r.origin_domain === "Bacteria");
  const otherRows = rows.filter(r => r.origin_domain !== "Bacteria");
  let html = "";
  if (bacteriaRows.length) {
    const rankOrder = ["phylum", "class_name", "order", "family", "genus", "species", "source_name"];
    const sortedBacteriaRows = bacteriaRows.slice().sort((a,b) => {
      const av = rankOrder.map(k => a[k] || "").join("|");
      const bv = rankOrder.map(k => b[k] || "").join("|");
      return av.localeCompare(bv, undefined, {numeric:true});
    });
    const body = sortedBacteriaRows.map(r => {
      const lineage = [
        r.phylum,
        r.class_name,
        r.order,
        r.family,
        r.genus,
        r.species || r.source_name
      ].filter(Boolean).join(" > ");
      return `
      <tr>
        <td>${esc(lineage || "lineage not specified")}</td>
        <td>${esc(r.ncbi_id || "")}</td>
      </tr>`;
    }).join("");
    html += `
      <h5>Bacteria lineages from MetOriginDB</h5>
      <p class="hint">Each row is one MetOrigin bacteria record, ordered from broad to specific: phylum > class > order > family > genus > species. Rows are not merged into a single taxonomic tree.</p>
      <table>
        <tr><th>Lineage</th><th>NCBI ID</th></tr>
        ${body}
      </table>`;
  }
  if (otherRows.length) {
    const root = {};
    for (const r of otherRows) {
      nestedInsert(root, [r.origin_domain || "Other origin", r.detail || r.source_name || "detail not specified"]);
    }
    html += `<h5>Other MetOrigin origin details</h5><div class="lineage-tree">${treeToHtml(root)}</div>`;
  }
  return html || "<p class='hint'>No expanded lineage rows available for this feature.</p>";
}
function renderLineageTreeOldMerged(rows) {
  const root = {};
  for (const r of rows) {
    if (r.origin_domain === "Bacteria") {
      nestedInsert(root, [
        "Bacteria",
        r.phylum || "phylum not specified",
        r.class_name || "class not specified",
        r.order || "order not specified",
        r.family || "family not specified",
        r.genus || "genus not specified",
        r.species || r.source_name || "species not specified"
      ]);
    } else {
      nestedInsert(root, [r.origin_domain || "Other origin", r.detail || r.source_name || "detail not specified"]);
    }
  }
  return `<div class="lineage-tree">${treeToHtml(root)}</div>`;
}
function renderSummary() {
  const rows = [
    ["Total features in report", DATA.summary.feature_count],
    ["MetOrigin matched identities", DATA.summary.metorigin_matched_features],
    ["domainMASST matched features", DATA.summary.domainmasst_matched_features],
    ["Integrated evidence features", DATA.summary.integrated_evidence_features],
    ["Features without any hit", DATA.summary.features_without_any_hit]
  ];
  let htmlRows = rows.map(r => `<tr><th>${esc(r[0])}</th><td>${esc(r[1])}</td></tr>`).join("");
  if (DATA.summary.metorigin_summary_rows.length) {
    htmlRows += `<tr><th colspan="2">MetOrigin matching summary</th></tr>`;
    for (const r of DATA.summary.metorigin_summary_rows) htmlRows += `<tr><td>${esc(Object.values(r).slice(0,2).join(" / "))}</td><td>${esc(Object.values(r).slice(-1)[0])}</td></tr>`;
  }
  document.getElementById("summaryTable").innerHTML = htmlRows;
}
function refresh() {
  const features = filteredFeatures();
  renderSankey(features);
  renderCircular(features);
}
function init() {
  document.getElementById("datasetName").textContent = DATA.meta.dataset;
  document.getElementById("analysisName").textContent = DATA.meta.analysis;
  document.getElementById("featureCount").textContent = DATA.summary.feature_count;
  document.getElementById("noHitCount").textContent = DATA.summary.features_without_any_hit;
  document.getElementById("metParams").innerHTML = tableFromObject(DATA.params.metorigin);
  document.getElementById("masstParams").innerHTML = tableFromObject(DATA.params.domainmasst);
  fillSelect("confidenceFilter", uniq(DATA.features.map(f => String(f.confidence_level))), "All confidence levels");
  fillSelect("evidenceFilter", uniq((DATA.evidence || []).map(e => e.evidence_source)), "MetOrigin + domainMASST");
  fillCheckboxGroup("originFilter", uniq((DATA.evidence || []).map(e => e.origin_domain)));
  fillSelect("sourceFilter", uniq(DATA.features.flatMap(f => [f.annotation_source, f.annotation_type])), "All annotation sources");
  fillSelect("classFilter", uniq(DATA.features.flatMap(f => [f.np_pathway, f.np_superclass, f.np_class])), "All NPClassifier classes");
  ["confidenceFilter","evidenceFilter","sourceFilter","classFilter","searchBox"].forEach(id => document.getElementById(id).addEventListener("input", refresh));
  document.getElementById("originFilter").addEventListener("change", refresh);
  renderSummary();
  refresh();
}
init();
</script>
</body>
</html>
"""
    out_path.write_text(template.replace("__TITLE__", title).replace("__DATA__", data_json), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", help="Run folder to report, e.g. runs/CURRENT_DATASET")
    parser.add_argument("--all-runs", action="store_true", help="Generate reports for every folder under runs/")
    args = parser.parse_args()

    project_dir = Path.cwd()
    if args.all_runs:
        run_dirs = [p for p in (project_dir / "runs").iterdir() if p.is_dir()]
    elif args.run:
        run_dirs = [Path(p) for p in args.run]
    else:
        raise SystemExit("Use --run runs/CURRENT_DATASET or --all-runs")

    for run_dir in run_dirs:
        run_dir = run_dir.resolve()
        if args.all_runs and not (run_dir / "outputs/metorigin_matches/bacterial_lineages.csv").exists():
            print(f"Keeping historical report unchanged: {run_dir.name}")
            continue
        payload = build_payload(run_dir, project_dir.resolve())
        out_path = run_dir / "interactive_origin_report.html"
        write_report(payload, out_path)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
