"""Read-only index validation summary; does not alter the original database CSV."""
import argparse
import json
from pathlib import Path
from metorigin_index import connect


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--index',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    db=connect(args.index,readonly=True)
    tables=['database_records','identifiers','source_associations','bacterial_lineages','lineage_definitions']
    result={name:db.execute('SELECT count(*) FROM '+name).fetchone()[0] for name in tables}
    result['issues']={r[0]:r[1] for r in db.execute('SELECT issue,count(*) FROM validation_issues GROUP BY issue')}
    result['identifier_validation']=[dict(r) for r in db.execute('SELECT namespace,validation,count(*) AS n FROM identifiers GROUP BY namespace,validation')]
    result['pubchem_cache']={r[0]:r[1] for r in db.execute('SELECT status,count(*) FROM pubchem_validation GROUP BY status')}
    result['metadata']={r[0]:json.loads(r[1]) for r in db.execute('SELECT * FROM metadata')}
    db.close()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='metadata'},indent=2))


if __name__=='__main__':main()
