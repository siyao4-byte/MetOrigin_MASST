import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from metorigin_index import build_index, connect, expand_lineages, candidates, structure_key, validate_cids
from indexed_match import match_row, run, input_scope_summary


class IndexedMetoriginTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / 'database.csv'
        self.index = self.root / 'index.sqlite'
        fields = ['Compound_name','PUBCHEM_COMPOUND_ID','PUBCHEM_SUBSTANCE_ID','SMILES_ID',
                  'from_human','from_which_part','from_bacteria','from_which_bacteria',
                  'bacteria_ncbi_id','bacteria_phylum','bacteria_class','bacteria_order',
                  'bacteria_family','bacteria_genus','bacteria_species','Database_source']
        self.rows = [dict(Compound_name='Example A',PUBCHEM_COMPOUND_ID='428{}4233',
            PUBCHEM_SUBSTANCE_ID='4233',SMILES_ID='NCCCN',from_human='Yes',from_which_part='Blood{}Urine',
            from_bacteria='Yes',from_which_bacteria='Species A{}Species B',bacteria_ncbi_id='1{}2',
            bacteria_phylum='Phylum P{}Phylum P',bacteria_family='Unknown{}Family B'),
            dict(Compound_name='Different compound',PUBCHEM_SUBSTANCE_ID='428',SMILES_ID='CCO'),
            dict(Compound_name='Example A alias',SMILES_ID='NCCCN',from_human='Yes',from_which_part='Feces')]
        with self.source.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(self.rows)
        build_index(self.source,self.index)
        self.db=connect(self.index)
        self.key=structure_key('NCCCN')
        validate_cids(self.db,['428','4233'],lambda ids:[{'CID':int(i),'InChIKey':self.key if i=='428' else structure_key('CCO')} for i in ids])

    def tearDown(self):
        self.db.close()
        self.tmp.cleanup()

    def test_synchronized_expansion(self):
        rows,issue=expand_lineages(self.rows[0])
        self.assertIsNone(issue)
        self.assertEqual([r[2] for r in rows],['Phylum P','Phylum P'])
        self.assertEqual([r[5] for r in rows],['','Family B'])

    def test_malformed_lengths_quarantined(self):
        rows,issue=expand_lineages(dict(self.rows[0],bacteria_phylum='Only one known phylum'))
        self.assertEqual(rows,[])
        self.assertEqual(issue['field'],'bacteria_phylum')

    def test_cid_never_queries_sid(self):
        self.assertNotIn(2,[r['record_id'] for r in candidates(self.db,'cid','428')])
        with self.assertRaises(ValueError): candidates(self.db,'sid','428')

    def test_wrong_database_cid_rejected(self):
        selected,_,evidence=match_row(self.db,{'cid':['4233']})
        self.assertEqual(selected,[])
        self.assertIn('database_cid_structure_conflict',[r['validation_status'] for r in evidence])

    def test_unvalidated_cid_is_not_sid_fallback(self):
        selected,_,_=match_row(self.db,{'cid':['999999']})
        self.assertEqual(selected,[])

    def test_compatible_lower_priority_evidence_retained(self):
        selected,method,_=match_row(self.db,{'cid':['428'],'smiles':['NCCCN'],'structure_key':self.key})
        self.assertEqual(selected,[1,3])
        self.assertEqual(method,'PubChem CID')

    def test_unknown_source_retained(self):
        row=self.db.execute("SELECT status FROM source_associations WHERE record_id=1 AND source='drug'").fetchone()
        self.assertEqual(row[0],'Unknown')

    def test_named_and_smiles_only_denominators(self):
        rows = [dict(compound_name='Example A', smiles='CCO', input_row_id='1', metorigin_matched=True),
                dict(compound_name='example a', smiles='CCO', input_row_id='2', metorigin_matched=False),
                dict(compound_name='NA', smiles='CCC', input_row_id='3', metorigin_matched=False),
                dict(compound_name='NA', smiles='NA', input_row_id='4', metorigin_matched=False)]
        scope = input_scope_summary(rows)
        self.assertEqual(scope['named'], dict(input_rows=2, unique_entries=1, matched_unique_entries=1))
        self.assertEqual(scope['smiles_only']['unique_entries'], 1)
        self.assertEqual(scope['unidentified']['unique_entries'], 1)

    def test_snapshot_and_all_associations(self):
        manifest=json.loads(self.index.with_suffix('.snapshot.json').read_text())
        self.assertTrue(manifest['complete'])
        self.assertEqual(self.db.execute('SELECT count(*) FROM bacterial_lineages WHERE record_id=1').fetchone()[0],2)
        self.assertEqual(self.db.execute("SELECT count(*) FROM source_associations WHERE record_id=1 AND source='human'").fetchone()[0],2)

    def test_end_to_end_exports_and_no_overwrite(self):
        (self.root/'config').mkdir()
        (self.root/'input.csv').write_text('compound.name,CID,smiles\nExample A,428,NCCCN\n',encoding='utf-8')
        config=self.root/'config/config.yaml'
        config.write_text('maps_annotation_file: input.csv\nmetorigin_index: index.sqlite\noutput_file: runs/test/outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv\n',encoding='utf-8')
        result=run(config,validate=False)
        summary=json.loads((result/'outputs/qc_reports/indexed_matching_summary.json').read_text())
        self.assertEqual(summary['matched_unique_identities'],1)
        self.assertEqual(summary['identity_source_links'],2)
        with (result/'outputs/metorigin_matches/bacterial_lineages.csv').open() as f:
            self.assertEqual(len(list(csv.DictReader(f))),2)
        with self.assertRaises(FileExistsError):run(config,validate=False)

    def test_display_limits_do_not_mutate_evidence(self):
        import generate_interactive_origin_report as report
        directory=self.root/'outputs/metorigin_matches'
        directory.mkdir(parents=True)
        fields=['feature_id','identity_id','compound_name','record_id','organism','ncbi_id','phylum','class_name','tax_order','family','genus','species']
        with (directory/'bacterial_lineages.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for i in range(50):w.writerow(dict(feature_id='f1',organism=f'Species {i}',phylum='P'))
        (directory/'source_associations.csv').write_text('feature_id,source,status,detail\nf1,human,Yes,Blood\nf1,human,Yes,Urine\n')
        original=(directory/'bacterial_lineages.csv').read_bytes()
        records=report.build_lineage_rows(self.root)
        self.assertEqual(sum(r['origin_domain']=='Bacteria' for r in records),41)
        self.assertEqual((directory/'bacterial_lineages.csv').read_bytes(),original)


if __name__ == '__main__':unittest.main()
