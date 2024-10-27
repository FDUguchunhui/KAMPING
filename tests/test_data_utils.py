from unittest.mock import Mock

import h5py
import numpy as np
import pytest

from kamping.parser import network, convert
from kamping.data.utils import convert_to_pyg, get_uniprot_protein_embeddings, convert_to_hetero_pyg, combine_graphs, \
    get_mol_embeddings, fetch_mol_file_string


def test_fetch_mol_file_string():
    mol_string = fetch_mol_file_string('C00001')
    # assert not none and not "Unhandled"
    assert mol_string != 'Unhandled'
    assert mol_string is not None


def test_get_mol_embedding():
    fake_instance = Mock()
    fake_instance.compounds = ['C00001', 'C00002']
    embeddings = get_mol_embeddings(fake_instance, transformer='morgan', radius=2)
    assert embeddings is not None


def test_get_protein_embedding(tmp_path):
    # simulate a h5 file
    embeddings = {'protein1': np.array([1, 2, 3]), 'protein2': np.array([4, 5, 6])}
    # save embeddings to pytest tmp
    with h5py.File(tmp_path / 'protein_embeddings.h5', 'w') as h5file:
        for key, value in embeddings.items():
            h5file.create_dataset(key, data=value)
    # test if the file is saved
    assert (tmp_path / 'protein_embeddings.h5').exists()
    # load the file
    # create mock object
    fake_instance = Mock()
    fake_instance.proteins = ['protein1', 'protein2']
    fake_instance.protein_id_type = 'uniprot'
    protein_embeddings = get_uniprot_protein_embeddings(fake_instance, tmp_path / 'protein_embeddings.h5')
    # assert protein_embeddings is dict of two
    assert protein_embeddings is not None
    assert len(protein_embeddings) == 2


class TestUtilsWorkWithGraph:
    converter = convert.Converter(species='hsa', gene_target='uniprot')

    def test_convert_to_pyg_gene_graph(self):
        gene_graph = network.KeggGraph('data/hsa00010.xml', type='gene',
                                       protein_group_as_interaction=True,
                                       multi_substrate_as_interaction=True)
        self.converter.convert(gene_graph)
        embeddings = get_uniprot_protein_embeddings(gene_graph, 'data/hsa00010_protein_embeddings.h5')
        data = convert_to_pyg(gene_graph, embeddings=embeddings)
        print(data)
        assert data is not None

    def test_convert_to_pyg_metablite_graph(self):
        metabolite_graph = network.KeggGraph('data/hsa00010.xml', type='metabolite',
                                             protein_group_as_interaction=True,
                                             multi_substrate_as_interaction=True)
        self.converter.convert(metabolite_graph)
        embeddings = {compound: np.zeros(1024) for compound in metabolite_graph.compounds}
        data = convert_to_pyg(metabolite_graph, embeddings=embeddings)
        print(data)
        assert data is not None

    def test_convert_to_hetero_pyg(self):
        mixed_graph = network.KeggGraph('data/hsa00010.xml', type='mixed',
                                        protein_group_as_interaction=True,
                                        multi_substrate_as_interaction=True)
        self.converter.convert(mixed_graph)

        protein_embeddings = {protein: np.zeros(1024) for protein in mixed_graph.proteins}
        mol_embeddings = {compound: np.zeros(1024) for compound in mixed_graph.compounds}
        data = convert_to_hetero_pyg(mixed_graph,
                                     protein_embeddings=protein_embeddings,
                                     mol_embeddings=mol_embeddings)
        print(data)
        assert data is not None


    def test_combine_graphs(self):
        graph1 = network.KeggGraph('data/hsa00010.xml', type='mixed',
                                        protein_group_as_interaction=True,
                                        multi_substrate_as_interaction=True)
        graph2 = network.KeggGraph('data/hsa00020.xml', type='mixed',
                                   protein_group_as_interaction=True,
                                   multi_substrate_as_interaction=True)
        self.converter.convert(graph1)
        self.converter.convert(graph2)
        protein_embeddings = {protein: np.zeros(1024) for protein in {*graph1.proteins, *graph2.proteins}}
        mol_embeddings = {compound: np.zeros(1024) for compound in {*graph1.compounds, *graph2.compounds}}
        data1 = convert_to_hetero_pyg(graph1,
                                     protein_embeddings=protein_embeddings,
                                     mol_embeddings=mol_embeddings)
        data2 = convert_to_hetero_pyg(graph2,
                                       protein_embeddings=protein_embeddings,
                                       mol_embeddings=mol_embeddings)
        print(data1)
        combined = combine_graphs([data1, data2])
        print(combined)
