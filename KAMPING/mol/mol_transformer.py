import contextlib
import io
import logging
from contextlib import contextmanager
import sys, os
import warnings
from typing import Literal, Any

import h5py
import numpy as np
warnings.filterwarnings("ignore", category=DeprecationWarning)
import scikit_mol.fingerprints
from rdkit.Chem import PandasTools
import pandas as pd


# supress warnings from RDKit DeprecationWarning and Parsing error
from rdkit import RDLogger
# Disable specific log
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')

class SmilesFileTransformer():
    '''
    Class to handle molecular vectorization
    '''

    def __init__(self, transformer: Any):
        self.transformer = transformer

    def transform(self, file:str, smiles_col:str, id_col:str) -> tuple[list[str], np.array]:
        '''
        Get molecular vector from MOL file with SMILES column and ID column

        params:
        file: str, path to the MOL file
        smiles_col: str, column name for SMILES
        id_col: str, column name for ID

        return:
        mol_vector: np.array, molecular vector

        usage:
        mol_vector = get_mol_vector('data/mol_files/smiles_file.tsv', smiles_col='smiles', id_col='id')
        '''

        df = pd.read_csv(file, sep='\t')

        with contextlib.redirect_stderr(None):
            PandasTools.AddMoleculeColumnToFrame(df, smilesCol=smiles_col)

        # get rows id with NaN in the ROMol column
        unvalid_row_id = df[df['ROMol'].isna()][id_col].tolist()

        logging.warning(f'Successfully parse {len(df) - len(unvalid_row_id)} rows with valid SMILES from the MOL file!\n'
                        f'total {len(unvalid_row_id)} Invalid rows with "Unhandled" in the ROMol column:\n'
                        f' {unvalid_row_id}, removed from the final output!')

        # remove rows with NaN in the ROMol column
        valid_row_id = df[~df['ROMol'].isna()][id_col].tolist()
        df = df.dropna(subset=['ROMol'])
        # get the molecular vector

        mol_embeddings = self.transformer.transform(df['ROMol'])
        return valid_row_id, mol_embeddings


class SmilesFileSMTransformer(SmilesFileTransformer):
    ''''
    A class to handle molecular vectorization using SMILES file using scikit-mol
    '''

    def __init__(self, transformer=Literal['morgan', 'rdkit', 'atom-pair', 'topological'], dim:int=2048,
                 **kwargs):
        if transformer == 'morgan':
            transformer = scikit_mol.fingerprints.MorganFingerprintTransformer(nBits=dim, **kwargs)
        elif transformer == 'rdkit':
            transformer = scikit_mol.fingerprints.RDKitFingerprintTransformer(fpSize=dim, **kwargs)
        elif transformer == 'atom-pair':
            transformer = scikit_mol.fingerprints.AtomPairFingerprintTransformer(nBits=dim, **kwargs)
        elif transformer == 'topological':
            transformer = scikit_mol.fingerprints.TopologicalTorsionFingerprintTransformer(nBits=dim, **kwargs)
        else:
            raise ValueError('Invalid transformer')

        super().__init__(transformer)

    def convert_mol_to_embed(self, file:str, smiles_col:str,
                             id_col:str,
                             output_path: str) -> None:
        '''
        Convert molecular files to SMILES

        params:
        file: str, path to the MOL file
        smiles_col: str, column name for SMILES
        id_col: str, column name for ID
        output_path: str, output path to save the h5 file

        return:
        None

        '''

        mol_ids, mol_embeddings = self.transform(file, smiles_col, id_col)
        # dictionary of key and embedding
        embeddings = dict(zip(mol_ids, mol_embeddings))

        # create parent directory if not exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # save  mol_ids, and mol_embedding as h5
        with h5py.File(output_path, 'w') as h5file:
            for key, value in embeddings.items():
                h5file.create_dataset(key, data=value)



# todo: another method to save embedding as pickle

if __name__ == '__main__':
    # supress warnings

    # set working directory use os
    # initialize the transformer
    mol_transformer = SmilesFileSMTransformer(transformer="rdkit", dim=2048)
    mol_transformer.convert_mol_to_embed( 'data/mol_files/smiles_file.tsv', 'smiles', 'id',
                                            output_path='data/embedding/mol_embedding.h5')


