import os

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
def read_all_files(directory):
    '''
    Read all files in the directory as one dataframe.
    '''
    # list a files except those start with "."

    all_files = os.listdir(directory)
    df = pd.DataFrame()
    for file in all_files:
        if file.startswith('.'):
            continue
        df = pd.concat([df, pd.read_csv(os.path.join(directory, file), sep='\t')])
    return df

# get unique values from column entry1 and entry2 combined and name starts with 'cpd'
def get_unique_compound_values(df):
    '''
    Get unique values from column entry1 and entry2 combined.
    '''
    all_entries = pd.concat([df['entry1'], df['entry2']]).unique()
    # convert all_entries to string

    compound_values = [entry for entry in all_entries if str(entry).startswith('cpd')]
    return compound_values
