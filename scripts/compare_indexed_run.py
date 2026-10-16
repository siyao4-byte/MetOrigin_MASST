"""Compare recorded associations on shared names without assuming identical input scope."""
import argparse
import csv
import json
import sys
from pathlib import Path
from metorigin_index import clean, normalize, sha256

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def read_identities(path):
    result = {}
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            name = row.get('compound_name',row.get('compound.name',''))
            key = normalize('name',name) or clean(row.get('smiles'))
            if not key:
                continue
            record = result.setdefault(key,{'matched':False,'sources':set()})
            record['matched'] |= str(row.get('metorigin_matched','')).lower() in ('true','yes','1')
            record['sources'].update(s for s in ('human','bacteria','plant','animal','food','drug','environment')
                                     if str(row.get('from_'+s,'')).lower() in ('yes','true','1'))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',required=True,type=Path)
    parser.add_argument('--current',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    old,new=read_identities(args.baseline),read_identities(args.current)
    shared=old.keys() & new.keys()
    rows=[]
    for key in sorted(shared):
        a,b=old[key],new[key]
        state='gained' if b['matched'] and not a['matched'] else 'lost' if a['matched'] and not b['matched'] else 'retained' if b['matched'] else 'unmatched'
        rows.append(dict(identity_key=key,change=state,old_sources=';'.join(sorted(a['sources'])),new_sources=';'.join(sorted(b['sources']))))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['identity_key','change','old_sources','new_sources']);w.writeheader();w.writerows(rows)
    summary={'baseline_sha256':sha256(args.baseline),'current_sha256':sha256(args.current),
             'baseline_identities':len(old),'current_identities':len(new),'shared_identities':len(shared),
             'counts':{s:sum(r['change']==s for r in rows) for s in ('gained','lost','retained','unmatched')},
             'source_sets_changed':sum(r['old_sources']!=r['new_sources'] for r in rows),
             'caveat':'Shared-name comparison only. Historical input checksum/version is not established; differences cannot be attributed solely to the matcher.'}
    args.output.with_suffix('.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
