"""Canonical indexed matching entry point. Total-list PubChem numbers are CID only."""
import re


def clean_colname(value):
    return re.sub(r'[^a-z0-9]+', '_', str(value).strip().lower()).strip('_')


def clean_dataframe_columns(df):
    result = df.copy()
    seen = {}
    columns = []
    for original in df.columns:
        name = clean_colname(original)
        count = seen.get(name, 0)
        columns.append(name if count == 0 else f'{name}_{count}')
        seen[name] = count + 1
    result.columns = columns
    return result


def detect_maps_columns(df):
    aliases = {
        'feature_id': ['feature_id', 'featureid', 'id'],
        'usi': ['feature_usi', 'usi'],
        'name': ['compound_name', 'name', 'title'],
        'smiles': ['smiles', 'canonical_smiles'],
        'inchikey': ['inchikey', 'inchi_key'],
        'cid': ['cid', 'pubchem_cid', 'pubchem_id', 'pubchem_compound_id'],
        'hmdb': ['hmdb_id', 'primary_hmdb_id', 'hmdb'],
        'kegg': ['kegg_id', 'kegg'],
    }
    return {key: next((name for name in names if name in df.columns), None)
            for key, names in aliases.items()}


if __name__ == '__main__':
    from indexed_match import main
    main()
