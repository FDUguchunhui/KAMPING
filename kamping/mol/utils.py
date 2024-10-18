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
