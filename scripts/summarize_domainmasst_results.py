#!/usr/bin/env python
"""
Summarize completed domainMASST batch outputs.

This is a lightweight reporting step for checking which MetOrigin-unmatched
features received domainMASST evidence, and which JSON tree subclasses are
represented among the matched nodes.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd


JSON_DOMAIN_SUFFIXES = ("combined", "microbe", "microbiome", "plant", "tissue", "food")


def infer_feature_id(path: Path) -> Optional[str]:
    match = re.search(r"AHGMD_unmatched_(.+?)_(?:combined|microbe|microbiome|plant|tissue|food|counts_|matches|datasets|library)", path.name)
    if match:
        return match.group(1)
    match = re.search(r"feature[_ -]?([A-Za-z0-9_.-]+)", path.name, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return None


def infer_json_domain(path: Path) -> str:
    stem = path.stem.lower()
    for domain in JSON_DOMAIN_SUFFIXES:
        if stem.endswith(f"_{domain}"):
            return domain
    return "unknown"


def infer_counts_domain(path: Path) -> str:
    match = re.search(r"_counts_([A-Za-z0-9]+)\.tsv$", path.name, flags=re.IGNORECASE)
    return match.group(1).lower() if match else "unknown"


def as_number(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def walk_json_tree(node: Any, path_parts: Optional[List[str]] = None) -> Iterable[Dict[str, Any]]:
    path_parts = path_parts or []
    if isinstance(node, dict):
        name = str(node.get("name") or node.get("label") or node.get("id") or "").strip()
        rank = str(node.get("Rank") or node.get("rank") or "").strip()
        current_path = path_parts + ([f"{rank}:{name}" if rank else name] if name else [])

        matched_size = as_number(node.get("matched_size"))
        if matched_size > 0 and name:
            yield {
                "subclass_name": name,
                "subclass_rank": rank if rank else "unranked",
                "matched_size": matched_size,
                "group_size": as_number(node.get("group_size")),
                "occurrence_fraction": as_number(node.get("occurrence_fraction")),
                "tree_path": " > ".join(part for part in current_path if part),
            }

        for child in node.get("children", []):
            yield from walk_json_tree(child, current_path)
    elif isinstance(node, list):
        for item in node:
            yield from walk_json_tree(item, path_parts)


def summarize_json_files(domainmasst_dir: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for json_path in sorted(domainmasst_dir.rglob("*.json")):
        feature_id = infer_feature_id(json_path)
        domain = infer_json_domain(json_path)
        try:
            data = json.loads(json_path.read_text(errors="replace"))
        except Exception as exc:
            rows.append({
                "feature_id": feature_id,
                "domain": domain,
                "json_file": str(json_path),
                "subclass_name": "JSON parse failed",
                "subclass_rank": "error",
                "matched_size": 0,
                "group_size": 0,
                "occurrence_fraction": 0,
                "tree_path": str(exc),
            })
            continue

        matched_nodes = list(walk_json_tree(data))
        if not matched_nodes:
            rows.append({
                "feature_id": feature_id,
                "domain": domain,
                "json_file": str(json_path),
                "subclass_name": "No matched JSON subclass",
                "subclass_rank": "none",
                "matched_size": 0,
                "group_size": 0,
                "occurrence_fraction": 0,
                "tree_path": "",
            })
            continue

        for item in matched_nodes:
            rows.append({
                "feature_id": feature_id,
                "domain": domain,
                "json_file": str(json_path),
                **item,
            })
    return pd.DataFrame(rows)


def summarize_counts_files(domainmasst_dir: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for tsv_path in sorted(domainmasst_dir.rglob("*_counts_*.tsv")):
        feature_id = infer_feature_id(tsv_path)
        domain = infer_counts_domain(tsv_path)
        try:
            n_rows = len(pd.read_csv(tsv_path, sep="\t"))
        except Exception:
            n_rows = 0
        rows.append({
            "feature_id": feature_id,
            "domain": domain,
            "counts_file": str(tsv_path),
            "counts_rows": n_rows,
        })
    return pd.DataFrame(rows)


def write_summaries(domainmasst_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    json_long = summarize_json_files(domainmasst_dir)
    counts_long = summarize_counts_files(domainmasst_dir)

    if json_long.empty:
        json_long = pd.DataFrame(columns=[
            "feature_id", "domain", "json_file", "subclass_name", "subclass_rank",
            "matched_size", "group_size", "occurrence_fraction", "tree_path",
        ])
    if counts_long.empty:
        counts_long = pd.DataFrame(columns=["feature_id", "domain", "counts_file", "counts_rows"])

    json_presence = (
        json_long.assign(has_json_match=lambda x: x["matched_size"].fillna(0).astype(float) > 0)
        .groupby(["feature_id", "domain"], dropna=False)
        .agg(
            json_files=("json_file", "nunique"),
            has_json_match=("has_json_match", "max"),
            json_matched_nodes=("has_json_match", "sum"),
            json_total_matched_size=("matched_size", "sum"),
        )
        .reset_index()
    )

    counts_presence = (
        counts_long.assign(has_counts_rows=lambda x: x["counts_rows"].fillna(0).astype(int) > 0)
        .groupby(["feature_id", "domain"], dropna=False)
        .agg(
            counts_files=("counts_file", "nunique"),
            has_counts_rows=("has_counts_rows", "max"),
            counts_rows=("counts_rows", "sum"),
        )
        .reset_index()
    )

    feature_domain = json_presence.merge(counts_presence, on=["feature_id", "domain"], how="outer")
    for col in ["json_files", "json_matched_nodes", "json_total_matched_size", "counts_files", "counts_rows"]:
        feature_domain[col] = feature_domain[col].fillna(0)
    for col in ["has_json_match", "has_counts_rows"]:
        feature_domain[col] = feature_domain[col].fillna(False).astype(bool)
    feature_domain["has_domain_evidence"] = feature_domain["has_json_match"] | feature_domain["has_counts_rows"]

    domain_counts = (
        feature_domain.groupby("domain", dropna=False)
        .agg(
            n_features_with_any_domain_output=("feature_id", "nunique"),
            n_features_with_json_match=("has_json_match", "sum"),
            n_features_with_counts_rows=("has_counts_rows", "sum"),
            n_features_with_domain_evidence=("has_domain_evidence", "sum"),
            total_json_matched_nodes=("json_matched_nodes", "sum"),
            total_counts_rows=("counts_rows", "sum"),
        )
        .reset_index()
        .sort_values(["n_features_with_domain_evidence", "n_features_with_any_domain_output", "domain"], ascending=[False, False, True])
    )

    subclass_top = (
        json_long[json_long["matched_size"].fillna(0).astype(float) > 0]
        .groupby(["domain", "subclass_rank", "subclass_name"], dropna=False)
        .agg(
            n_features=("feature_id", "nunique"),
            total_matched_size=("matched_size", "sum"),
            median_occurrence_fraction=("occurrence_fraction", "median"),
            max_occurrence_fraction=("occurrence_fraction", "max"),
            example_feature_ids=("feature_id", lambda x: " | ".join(sorted(set(x.dropna().astype(str)))[:10])),
        )
        .reset_index()
        .sort_values(["domain", "n_features", "total_matched_size"], ascending=[True, False, False])
    )

    broad_taxonomy_ranks = {"domain", "superkingdom", "kingdom", "phylum", "class"}
    generic_names = {"root", "unknown", "origin", "host", "health phenotype"}
    subclass_top_curated = (
        subclass_top[
            (subclass_top["domain"] != "combined")
            & (~subclass_top["subclass_rank"].str.lower().isin(broad_taxonomy_ranks))
            & (~subclass_top["subclass_name"].str.lower().isin(generic_names))
        ]
        .sort_values(["domain", "n_features", "total_matched_size"], ascending=[True, False, False])
    )

    json_long.to_csv(out_dir / "domainmasst_testing_summary_json_subclasses_long.csv", index=False)
    counts_long.to_csv(out_dir / "domainmasst_testing_summary_counts_files_long.csv", index=False)
    feature_domain.to_csv(out_dir / "domainmasst_testing_summary_feature_domain_presence.csv", index=False)
    domain_counts.to_csv(out_dir / "domainmasst_testing_summary_domain_counts.csv", index=False)
    subclass_top.to_csv(out_dir / "domainmasst_testing_summary_json_subclasses_top.csv", index=False)
    subclass_top_curated.to_csv(out_dir / "domainmasst_testing_summary_json_subclasses_top_curated.csv", index=False)

    print("Domain counts")
    print(domain_counts.to_string(index=False))
    print()
    print("Top JSON subclasses")
    print(subclass_top_curated.head(30).to_string(index=False))
    print()
    print(f"Wrote summary CSVs to: {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize domainMASST JSON and counts outputs.")
    parser.add_argument("--domainmasst_dir", default="outputs/domainmasst_unmatched_ms2/domainMASST")
    parser.add_argument("--out_dir", default="outputs/domainmasst_unmatched_ms2/testing_summary")
    args = parser.parse_args()
    write_summaries(Path(args.domainmasst_dir), Path(args.out_dir))


if __name__ == "__main__":
    main()
