"""Add explicit named/SMILES-only denominators to an existing indexed run audit."""
import argparse
import csv
import json
from pathlib import Path
from indexed_match import input_scope_summary

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', required=True, type=Path)
    args = parser.parse_args()
    source = args.run / 'outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv'
    with source.open(encoding='utf-8', newline='') as handle:
        scope = input_scope_summary(csv.DictReader(handle))
    summary = args.run / 'outputs/qc_reports/indexed_matching_summary.json'
    result = json.loads(summary.read_text(encoding='utf-8'))
    result['input_scope_by_identity_type'] = scope
    summary.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(scope, indent=2))
