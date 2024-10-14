import pandas as pd
import requests
from rdkit import Chem

def fetch_mol_file_string(compound_id):
    '''
    Fetch the MOL file string for a given compound ID.
    '''
    urlStr = "https://www.genome.jp/dbget-bin/www_bget?-f+m+compound+" + compound_id
    molFile = requests.get(urlStr)
    return(str(molFile.text))

def get_smiles(mol_file_string):
    # Use RDKit to convert the mol file string to a SMILES.
    mol = Chem.MolFromMolBlock(mol_file_string);
    try:
        smiles_string = Chem.rdmolfiles.MolToSmiles(mol);
    except:
        return("Unhandled");        # This was a polymer, which RDKit can't handle.
    return(smiles_string);


# read all file in the directory as one dataframe

# get unique values from column entry1 and entry2 combined and name starts with 'cpd'
def get_unique_compound_values(df):
    '''
    Get unique values from column entry1 and entry2 combined.
    '''
    # get all emtries from entry1 if entry1_type is compound
    # get all entries from entry2 if entry2_type is compound
    entry1_compound = df[df['entry1_type'] == 'compound']['entry1']
    entry2_compound = df[df['entry2_type'] == 'compound']['entry2']
    compounds = pd.concat([entry1_compound, entry2_compound]).unique()

    return compounds
