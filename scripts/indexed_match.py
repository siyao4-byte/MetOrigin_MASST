"""Indexed MAPS matching with complete, uncapped association exports."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import yaml
from metorigin_index import SOURCES, candidates, clean, connect, normalize, sha256, structure_key, tokens, validate_cids

LABELS = {'cid': 'PubChem CID', 'smiles': 'SMILES exact', 'inchikey': 'InChIKey',
          'hmdb': 'HMDB ID', 'kegg': 'KEGG ID', 'name': 'compound name exact'}


def match_row(db, keys):
    query_structure = keys.get('structure_key', '')
    evidence = []
    accepted = {}
    for priority, namespace in enumerate(('cid', 'smiles', 'inchikey', 'hmdb', 'kegg', 'name')):
        for value in keys.get(namespace, []):
            for record in candidates(db, namespace, value):
                rid, rkey = record['record_id'], record['structure_key']
                state, ok = 'source_identifier_match', True
                if query_structure and rkey and query_structure != rkey:
                    state, ok = 'structure_conflict', False
                elif query_structure and rkey:
                    state = 'exact_structure_agreement'
                if namespace == 'cid' and ok:
                    validated = db.execute('SELECT inchikey,status FROM pubchem_validation WHERE cid=?', (value,)).fetchone()
                    if validated and validated['status'] == 'resolved' and rkey:
                        if validated['inchikey'] != rkey:
                            state, ok = 'database_cid_structure_conflict', False
                        elif query_structure and validated['inchikey'] != query_structure:
                            state, ok = 'input_cid_structure_conflict', False
                        else:
                            state = 'pubchem_structure_validated'
                    elif not (query_structure and rkey and query_structure == rkey):
                        state, ok = 'cid_unverified_not_used', False
                if namespace == 'name' and ok and state == 'source_identifier_match':
                    state = 'name_only_unverified'
                evidence.append({'record_id': rid, 'namespace': namespace, 'value': value,
                                 'validation_status': state, 'compatible': ok, 'selected': False})
                if ok:
                    accepted[rid] = min(priority, accepted.get(rid, priority))
    if not accepted:
        return [], 'unmatched', evidence
    best = min(accepted.values())
    anchors = [rid for rid, rank in accepted.items() if rank == best]
    anchor_keys = {db.execute('SELECT structure_key FROM database_records WHERE record_id=?', (rid,)).fetchone()[0] for rid in anchors}
    anchor_keys.discard('')
    # Include lower-priority records only when they agree with an anchored exact structure.
    selected = set(anchors)
    for rid in accepted:
        rkey = db.execute('SELECT structure_key FROM database_records WHERE record_id=?', (rid,)).fetchone()[0]
        if len(anchor_keys) == 1 and rkey in anchor_keys:
            selected.add(rid)
    if len(anchor_keys) > 1 and not query_structure:
        for e in evidence:
            if e['record_id'] in selected and e['compatible']:
                e['validation_status'] = 'ambiguous_multiple_structures'
        return [], 'ambiguous', evidence
    for e in evidence:
        e['selected'] = e['compatible'] and e['record_id'] in selected
    return sorted(selected), list(LABELS.values())[best], evidence


def write_csv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def input_scope_summary(rows):
    groups = {'named': set(), 'smiles_only': set(), 'unidentified': set()}
    matched = {key: set() for key in groups}
    counts = Counter()
    for row in rows:
        name = normalize('name', row.get('compound_name', ''))
        smiles = clean(row.get('smiles', ''))
        kind = 'named' if name else 'smiles_only' if smiles else 'unidentified'
        key = name or smiles or row['input_row_id']
        groups[kind].add(key)
        counts[kind] += 1
        if str(row.get('metorigin_matched', '')).lower() in ('true', '1', 'yes'):
            matched[kind].add(key)
    return {kind: {'input_rows': counts[kind], 'unique_entries': len(groups[kind]),
                   'matched_unique_entries': len(matched[kind])} for kind in groups}


def run(config_path, validate=True):
    from metorigin_match import clean_dataframe_columns, detect_maps_columns
    config_path = Path(config_path).resolve()
    project = config_path.parent.parent
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    if config.get('input_pubchem_namespace', 'cid') != 'cid' or config.get('allow_sid_fallback', False):
        raise ValueError('Total-list PubChem identifiers must be CID; SID fallback is forbidden.')
    resolve = lambda value: (project / value).resolve()
    input_path = resolve(config.get('maps_annotation_file', config.get('maps_file')))
    index_path = resolve(config['metorigin_index'])
    output = resolve(config['output_file'])
    outdir = output.parent
    if output.exists():
        raise FileExistsError('Run output already exists. Select a new run directory; historical results are not overwritten.')
    db = connect(index_path)
    meta = {r['key']: json.loads(r['value']) for r in db.execute('SELECT * FROM metadata')}
    if not meta.get('complete'):
        raise ValueError('Incomplete database index')
    snapshot = Path(meta['source_path'])
    if sha256(snapshot) != meta['source_sha256']:
        raise ValueError('Database CSV no longer matches the indexed snapshot')
    data = clean_dataframe_columns(pd.read_csv(input_path, dtype=str, keep_default_na=False))
    cols = detect_maps_columns(data)
    rows = data.to_dict('records')
    feature_values = [row.get(cols['feature_id'], '') for row in rows] if cols['feature_id'] else []
    preserve_feature_ids = bool(feature_values) and all(feature_values) and len(set(feature_values)) == len(feature_values)
    cids = set()
    keys_all = []
    for row in rows:
        keys = {}
        for namespace in LABELS:
            column = cols.get(namespace)
            raw = row.get(column, '') if column else ''
            values = [raw] if namespace in ('name', 'smiles') else tokens(raw)
            keys[namespace] = [normalize(namespace, v) for v in values if normalize(namespace, v)]
        skey = structure_key(keys['smiles'][0]) if keys['smiles'] else ''
        if not skey and len(keys['inchikey']) == 1:
            skey = keys['inchikey'][0]
        keys['structure_key'] = skey
        if skey:
            keys['inchikey'] = list(dict.fromkeys(keys['inchikey'] + [skey]))
        cids.update(keys['cid'])
        keys_all.append(keys)
    if validate:
        validate_cids(db, cids)
    outdir.mkdir(parents=True, exist_ok=True)
    run_dir = outdir.parent.parent if outdir.parent.name == 'outputs' else outdir
    qc = run_dir / 'outputs/qc_reports'
    qc.mkdir(parents=True, exist_ok=True)
    (run_dir / 'config').mkdir(exist_ok=True)
    shutil.copy2(config_path, run_dir / 'config' / config_path.name)
    edb = connect(outdir / 'match_evidence.sqlite')
    edb.execute('CREATE TABLE match_evidence(input_row_id TEXT, identity_id TEXT, record_id INTEGER, namespace TEXT, value TEXT, validation_status TEXT, compatible INTEGER, selected INTEGER)')
    handles = {}
    fields = {
        'source_associations': ['input_row_id','feature_id','identity_id','compound_name','record_id','source','status','detail','provenance'],
        'bacterial_lineages': ['input_row_id','feature_id','identity_id','compound_name','record_id','organism','ncbi_id','phylum','class_name','tax_order','family','genus','species'],
        'match_evidence': ['input_row_id','identity_id','record_id','namespace','value','validation_status','compatible','selected'],
    }
    for name, columns in fields.items():
        handle = (outdir / (name + '.csv')).open('w', encoding='utf-8', newline='')
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore'); writer.writeheader()
        handles[name] = (handle, writer)
    matched_ids, features_by_source, associations = set(), defaultdict(set), set()
    exported_records = set()
    match_cache = {}
    record_status_cache = {}
    method_counts, validation_counts, conflicts = Counter(), Counter(), 0
    final = []
    for position, (row, keys) in enumerate(zip(rows, keys_all), 1):
        input_id = f'row_{position:07d}'
        name = row.get(cols['name'], '') if cols['name'] else ''
        identity_key = normalize('name', name) or (keys['smiles'][0] if keys['smiles'] else input_id)
        identity_id = 'identity_' + hashlib.sha256(identity_key.encode()).hexdigest()[:24]
        feature_id = row.get(cols['feature_id'], '') if cols['feature_id'] else ''
        row['original_feature_id'] = feature_id
        # Total lists can reuse feature IDs across assays. Row IDs are globally unique within the input.
        public_feature_id = feature_id if preserve_feature_ids else input_id
        row.update(input_row_id=input_id, feature_id=public_feature_id, identity_id=identity_id, compound_name=name)
        cache_key = json.dumps(keys, sort_keys=True)
        if cache_key not in match_cache:
            match_cache[cache_key] = match_row(db, keys)
        selected, method, evidence = match_cache[cache_key]
        for e in evidence:
            entry = dict(input_row_id=input_id, identity_id=identity_id, **e)
            handles['match_evidence'][1].writerow(entry)
            edb.execute('INSERT INTO match_evidence VALUES(?,?,?,?,?,?,?,?)', tuple(entry[k] for k in fields['match_evidence']))
            validation_counts[e['validation_status']] += 1
        conflicts += any('conflict' in e['validation_status'] for e in evidence)
        flags = defaultdict(set)
        for rid in selected:
            export_this = (identity_id, rid) not in exported_records
            base = dict(input_row_id='', feature_id=identity_id, identity_id=identity_id, compound_name=name, record_id=rid)
            if not export_this:
                for source, states in record_status_cache[rid].items():
                    flags[source].update(states)
                continue
            record_flags = defaultdict(set)
            for assoc in db.execute('SELECT * FROM source_associations WHERE record_id=?', (rid,)):
                a = dict(assoc)
                handles['source_associations'][1].writerow(dict(base, **{k:v for k,v in a.items() if k != 'record_id'}))
                flags[a['source']].add(a['status'])
                record_flags[a['source']].add(a['status'])
            record_status_cache[rid] = record_flags
            for lineage in db.execute('''SELECT d.* FROM bacterial_lineages b JOIN lineage_definitions d USING(lineage_id)
                                        WHERE b.record_id=? AND EXISTS(SELECT 1 FROM source_associations a
                                        WHERE a.record_id=b.record_id AND a.source='bacteria' AND a.status='Yes')''', (rid,)):
                handles['bacterial_lineages'][1].writerow(dict(base, **dict(lineage)))
            exported_records.add((identity_id, rid))
        for source in SOURCES:
            status = 'Yes' if 'Yes' in flags[source] else 'Unknown' if not flags[source] or 'Unknown' in flags[source] else 'No'
            row['from_' + source] = status
        positive = [s for s in SOURCES if row['from_' + s] == 'Yes']
        for source in positive:
            features_by_source[source].add(identity_id)
            associations.add((identity_id,source))
        row.update(metorigin_matched=bool(selected), metorigin_match_method=method,
                   metorigin_n_database_records=len(selected), metorigin_origin_categories=';'.join(positive),
                   metorigin_n_origin_categories=len(positive), human_only_recorded_association=positive == ['human'],
                   association_conflict=any(v >= {'Yes','No'} for v in flags.values()))
        # Details live in normalized tables, not repeated wide strings on every assay row.
        if selected:
            matched_ids.add(identity_id)
        method_counts[method] += 1
        final.append(row)
        if position % 1000 == 0:
            edb.commit()
            print(f'Matched {position}/{len(rows)} input rows', flush=True)
    edb.execute('CREATE INDEX match_input ON match_evidence(input_row_id,selected)')
    edb.execute('CREATE INDEX match_identity ON match_evidence(identity_id,selected)')
    edb.commit(); edb.close()
    for handle, _ in handles.values():
        handle.close()
    output_fields = list(dict.fromkeys(k for r in final for k in r))
    write_csv(output, final, output_fields)
    for suffix in ('_FULL', '_COMPACT'):
        shutil.copy2(output, output.with_name(output.stem + suffix + '.csv'))
    write_csv(qc / 'metorigin_unmatched_features.csv', [r for r in final if not r['metorigin_matched']], output_fields)
    summary = {'input_rows': len(final), 'unique_input_identities': len({r['identity_id'] for r in final}),
               'input_scope_by_identity_type': input_scope_summary(final),
               'matched_input_rows': sum(r['metorigin_matched'] for r in final),
               'matched_unique_identities': len(matched_ids), 'identity_source_links': len(associations),
               'unique_identities_by_source': {s:len(v) for s,v in features_by_source.items()},
               'match_methods_input_rows': dict(method_counts), 'candidate_validation': dict(validation_counts),
               'input_rows_with_rejected_conflicts': conflicts,
               'input_sha256': sha256(input_path), 'input_path': str(input_path),
               'original_feature_ids_preserved': preserve_feature_ids,
               'association_export_unit': 'identity_id x database record; all input links retained in match_evidence',
               'index_snapshot': meta, 'input_pubchem_namespace': 'CID only; no SID lookup or fallback',
               'interpretation': 'Recorded source associations, not definitive biological origins.'}
    (qc / 'indexed_matching_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(qc / 'metorigin_matching_summary.csv', [{'metorigin_match_method':m,'n_input_rows':n} for m,n in method_counts.items()], ['metorigin_match_method','n_input_rows'])
    db.close()
    print(json.dumps({k:v for k,v in summary.items() if k != 'index_snapshot'}, indent=2), flush=True)
    return run_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--offline', action='store_true', help='Use cached CID validation only; never reinterpret a CID as SID.')
    args = parser.parse_args()
    run(args.config, not args.offline)


if __name__ == '__main__':
    main()
