import os

import pytest
import KAMPING.mol.mol_transformer as mol_transformer

class TestMolTransformer:

    def test_smiles_file_transformer(self):
        # print current working directory
        print('Current working directory:', os.getcwd())

        transformer = mol_transformer.SmilesFileSMTransformer(transformer='morgan')
        compound_id,  embedding = transformer.transform('data/smiles_file.tsv', smiles_col='smiles', id_col='id')
        assert compound_id is not None
        assert embedding is not None
        # check len of compound_id and embedding first dimension is the same
        assert len(compound_id) == len(embedding)


    def test_smiles_file_transformer_with_invalid_transformer(self):
        with pytest.raises(ValueError, match='Invalid transformer'):
            mol_transformer.SmilesFileSMTransformer(transformer='invalid_transformer')

    def test_smile_file_transformer_snapshot(self, snapshot):
        transformer = mol_transformer.SmilesFileSMTransformer(transformer='morgan')
        compound_id, embedding = transformer.transform('data/smiles_file.tsv', smiles_col='smiles', id_col='id')
        snapshot.assert_match(embedding)