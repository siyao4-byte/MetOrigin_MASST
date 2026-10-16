"""Lossless, snapshot-addressed MetOriginDB index and typed matching primitives."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import time
import urllib.error
import urllib.request
import zlib
from functools import lru_cache
from pathlib import Path

SCHEMA_VERSION = 1
SOURCES = ('human', 'bacteria', 'plant', 'animal', 'environment', 'drug', 'food')
RANKS = ('bacteria_ncbi_id', 'bacteria_phylum', 'bacteria_class', 'bacteria_order',
         'bacteria_family', 'bacteria_genus', 'bacteria_species')
NULLS = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', '-1', '0'}


def clean(value):
    s = str(value or '').strip()
    return '' if s.casefold() in NULLS else s


def tokens(value):
    return list(dict.fromkeys(clean(x) for x in re.split(r'\{\}|\s*\|\|?\s*|;', str(value or '')) if clean(x)))


def normalize(namespace, value):
    value = clean(value)
    if namespace in ('cid', 'sid'):
        value = re.sub(r'\.0$', '', value)
        return value if re.fullmatch(r'[1-9][0-9]*', value) else ''
    if namespace == 'name':
        return ' '.join(value.casefold().split())
    if namespace in ('inchikey', 'hmdb', 'kegg'):
        return value.upper()
    return value


@lru_cache(maxsize=100000)
def structure_key(smiles):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(mol) if mol is not None else ''
    except Exception:
        return ''


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def connect(path, readonly=False):
    path = Path(path).resolve()
    db = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True) if readonly else sqlite3.connect(path)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA foreign_keys=ON')
    return db


def expand_lineages(row):
    # Keep positions and Unknown placeholders until whole lineage rows exist.
    names = str(row.get('from_which_bacteria') or '').split('{}')
    if not any(clean(x) for x in names):
        return [], None
    columns = [names]
    for field in RANKS:
        value = row.get(field, '')
        parts = str(value or '').split('{}') if clean(value) else [''] * len(names)
        if len(parts) != len(names):
            return [], {'field': field, 'organisms': len(names), 'values': len(parts)}
        columns.append(parts)
    return list(dict.fromkeys(tuple(clean(x) for x in lineage) for lineage in zip(*columns))), None


def build_index(source, target):
    source, target = Path(source).resolve(), Path(target).resolve()
    if target.exists():
        raise FileExistsError(f'Index exists; use a new snapshot path: {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    before = source.stat()
    checksum = sha256(source)
    db = connect(target)
    db.execute('PRAGMA journal_mode=WAL')
    db.executescript('''
      CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
      CREATE TABLE database_records(record_id INTEGER PRIMARY KEY, compound_name TEXT,
        structure_key TEXT, raw_json_zlib BLOB NOT NULL);
      CREATE TABLE identifiers(record_id INTEGER REFERENCES database_records, namespace TEXT,
        value TEXT, field TEXT, validation TEXT,
        PRIMARY KEY(record_id,namespace,value,field)) WITHOUT ROWID;
      CREATE TABLE source_associations(record_id INTEGER REFERENCES database_records,
        source TEXT, status TEXT, detail TEXT, provenance TEXT,
        PRIMARY KEY(record_id,source,status,detail)) WITHOUT ROWID;
      CREATE TABLE lineage_definitions(lineage_id INTEGER PRIMARY KEY, organism TEXT,
        ncbi_id TEXT, phylum TEXT, class_name TEXT, tax_order TEXT, family TEXT, genus TEXT, species TEXT,
        UNIQUE(organism,ncbi_id,phylum,class_name,tax_order,family,genus,species));
      CREATE TABLE bacterial_lineages(record_id INTEGER REFERENCES database_records,
        lineage_id INTEGER REFERENCES lineage_definitions, PRIMARY KEY(record_id,lineage_id)) WITHOUT ROWID;
      CREATE TABLE validation_issues(record_id INTEGER, issue TEXT, detail TEXT);
      CREATE TABLE pubchem_validation(cid TEXT PRIMARY KEY, inchikey TEXT, status TEXT, checked_utc TEXT);
    ''')
    meta = {'schema_version': SCHEMA_VERSION, 'source_path': str(source), 'source_sha256': checksum,
            'source_bytes': before.st_size, 'source_mtime_ns': before.st_mtime_ns,
            'builder_sha256': sha256(__file__), 'upstream': 'https://github.com/tidymass/metorigindb',
            'upstream_revision': 'not recorded in downloaded CSV', 'complete': False}
    fields = {'cid': ['PUBCHEM_COMPOUND_ID'], 'sid': ['PUBCHEM_SUBSTANCE_ID'],
              'inchikey': ['INCHIKEY_ID', 'INCHIKEY_ID_all'], 'hmdb': ['HMDB_ID', 'HMDB_ID_all'],
              'kegg': ['KEGG_ID', 'KEGG_ID_all'], 'smiles': ['SMILES_ID', 'SMILES_ID_all'],
              'name': ['Compound_name', 'Synonyms']}
    lineages_seen = {}
    csv.field_size_limit(100000000)
    n = 0
    with source.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        required = {'Compound_name', 'PUBCHEM_COMPOUND_ID', 'from_which_bacteria'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Unexpected MetOriginDB schema')
        meta['columns'] = reader.fieldnames
        for n, row in enumerate(reader, 1):
            skey = structure_key(clean(row.get('SMILES_ID'))) if clean(row.get('SMILES_ID')) else ''
            db.execute('INSERT INTO database_records VALUES(?,?,?,?)',
                       (n, row.get('Compound_name'), skey, zlib.compress(json.dumps(row).encode(), 1)))
            for ns, names in fields.items():
                for field in names:
                    for value in tokens(row.get(field)):
                        key = normalize(ns, value)
                        if not key:
                            db.execute('INSERT INTO validation_issues VALUES(?,?,?)', (n, 'invalid_identifier', field + ':' + value))
                            continue
                        status = 'namespace_declared_unverified' if ns in ('cid', 'sid') else 'source_declared'
                        db.execute('INSERT OR IGNORE INTO identifiers VALUES(?,?,?,?,?)', (n, ns, key, field, status))
            if skey:
                db.execute('INSERT OR IGNORE INTO identifiers VALUES(?,?,?,?,?)',
                           (n, 'inchikey', skey, 'SMILES_ID:RDKit', 'structure_derived'))
            for source_name in SOURCES:
                flag = clean(row.get('from_' + source_name)).casefold()
                status = 'Yes' if flag in ('yes', 'true', '1') else 'No' if flag in ('no', 'false') else 'Unknown'
                field = 'from_which_part' if source_name == 'human' else 'from_which_' + source_name
                for detail in tokens(row.get(field)) or ['']:
                    db.execute('INSERT OR IGNORE INTO source_associations VALUES(?,?,?,?,?)',
                               (n, source_name, status, detail, row.get('Database_source', '')))
            lineages, issue = expand_lineages(row)
            if issue:
                db.execute('INSERT INTO validation_issues VALUES(?,?,?)', (n, 'taxonomy_length_mismatch', json.dumps(issue)))
            for lineage in lineages:
                lid = lineages_seen.get(lineage)
                if lid is None:
                    lid = len(lineages_seen) + 1
                    db.execute('INSERT INTO lineage_definitions VALUES(?,?,?,?,?,?,?,?,?)', (lid,) + lineage)
                    lineages_seen[lineage] = lid
                db.execute('INSERT OR IGNORE INTO bacterial_lineages VALUES(?,?)', (n, lid))
            if n % 10000 == 0:
                db.commit()
                print(f'Indexed {n} records; {len(lineages_seen)} distinct complete lineages', flush=True)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('Source changed while indexing; index is not complete')
    db.executescript('''
      CREATE INDEX identifier_lookup ON identifiers(namespace,value);
      CREATE INDEX association_source ON source_associations(source,status);
      CREATE INDEX bacterial_taxon ON bacterial_lineages(lineage_id);
    ''')
    meta.update(complete=True, records=n)
    db.executemany('INSERT INTO metadata VALUES(?,?)', [(k, json.dumps(v)) for k, v in meta.items()])
    db.commit()
    assert db.execute('PRAGMA quick_check').fetchone()[0] == 'ok'
    db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
    db.close()
    target.with_suffix('.snapshot.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print(json.dumps(meta, indent=2), flush=True)


def validate_cids(db, cids, fetch=None):
    """Look up CID properties only. Never call a SID endpoint for input numbers."""
    def fetch_default(batch):
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + ','.join(batch) + '/property/InChIKey/JSON'
        with urllib.request.urlopen(url, timeout=45) as response:
            return json.load(response)['PropertyTable']['Properties']
    fetch = fetch or fetch_default
    pending = [x for x in sorted(set(cids)) if normalize('cid', x) and not db.execute(
        "SELECT 1 FROM pubchem_validation WHERE cid=? AND status IN ('resolved','not_found')", (x,)).fetchone()]
    def request(batch):
        try:
            rows = fetch(batch)
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404) and len(batch) > 1:
                middle = len(batch) // 2
                request(batch[:middle]); request(batch[middle:])
                return
            state = 'not_found' if exc.code == 404 else 'unavailable'
            rows = []
        except Exception:
            state, rows = 'unavailable', []
        else:
            state = 'not_found'
        found = {str(r['CID']): r.get('InChIKey', '') for r in rows}
        for cid in batch:
            db.execute('INSERT OR REPLACE INTO pubchem_validation VALUES(?,?,?,?)',
                       (cid, found.get(cid, ''), 'resolved' if cid in found else state, time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
            if cid in found:
                db.execute('''UPDATE identifiers SET validation=CASE
                    WHEN (SELECT structure_key FROM database_records r WHERE r.record_id=identifiers.record_id)=?
                      THEN 'pubchem_structure_validated'
                    WHEN (SELECT structure_key FROM database_records r WHERE r.record_id=identifiers.record_id)!=''
                      THEN 'pubchem_structure_conflict' ELSE 'cid_exists_structure_unchecked' END
                    WHERE namespace='cid' AND value=?''', (found[cid],cid))
        db.commit()
    for offset in range(0, len(pending), 100):
        request(pending[offset:offset+100])
        print(f'CID validation {min(offset+100,len(pending))}/{len(pending)}', flush=True)
        time.sleep(.3)


def candidates(db, namespace, value):
    # SID is intentionally absent from the supported input namespaces.
    if namespace not in ('cid', 'inchikey', 'smiles', 'hmdb', 'kegg', 'name'):
        raise ValueError('Unsupported input namespace; total-list PubChem values are CID only')
    return db.execute('''SELECT DISTINCT r.record_id,r.compound_name,r.structure_key
        FROM identifiers i JOIN database_records r USING(record_id)
        WHERE i.namespace=? AND i.value=?''', (namespace, normalize(namespace, value))).fetchall()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--index', required=True, type=Path)
    args = parser.parse_args()
    build_index(args.source, args.index)
