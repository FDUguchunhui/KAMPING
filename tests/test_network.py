import io
from enum import Enum
from unittest.mock import patch
import pytest
import pandas as pd
from pathlib import Path

import requests

from kamping.parser.convert import Converter
from kamping.parser.network import KeggGraph
from kamping.main import network
import xml.etree.ElementTree as ET
from click.testing import CliRunner
from unittest.mock import Mock
import networkx as nx


def parse_kgml_file(file_path, **kwargs):
    interaction = KeggGraph(input_data=file_path, **kwargs)
    return interaction.interaction


test_file = 'data/hsa00010_test.xml'


class TestGenesInteractionParser:
    xml_data = '''
    <root>
        <entry id="1" name="gene1" type="gene"/>
        <entry id="2" name="gene2" type="gene"/>
        <entry id="3" name="gene3" type="compound"/>
        <relation entry1="1" entry2="2" type="PPrel">
            <subtype name="compound" value="3"/>
        </relation>
    </root>
    '''
    # xml_data to file-like object
    xml_data = io.StringIO(xml_data)

    def test_parse_empty_kgml_file(self):
        with pytest.raises(OSError):
            parse_kgml_file('data/kgml_hsa/empty.xml', type='gene')

    def test_genes_parser_single_file(self):
        output_df = parse_kgml_file(test_file, type="gene", unique=False, verbose=False)
        assert not output_df.empty
        assert 'entry1' in output_df.columns
        assert 'entry2' in output_df.columns

    def test_genes_parser_directory(self):

        output_dir = Path('output')
        output_dir.mkdir(parents=True, exist_ok=True)
        runner = CliRunner()
        runner.invoke(network, ['--type', 'gene', 'hsa', test_file, '--out_dir', 'output'])
        output_files = list(output_dir.glob('*.tsv'))
        assert len(output_files) > 0
        for file in output_files:
            df = pd.read_csv(file, sep='\t')
            assert not df.empty
            assert 'entry1' in df.columns
            assert 'entry2' in df.columns

    def test_parse_without_id_conversion(self, snapshot):
        ''''
        Test without id_conversion
        '''
        output_df = parse_kgml_file(test_file, type="mpi", unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_with_id_conversion(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        id_converter = Converter('hsa', target='uniprot')
        output_df = parse_kgml_file(test_file, type="gene", id_converter=id_converter,
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_with_multiple_subtype(self):
        '''expect multiple rows for each relation if there are multiple subtypes'''
        xml_string = '''
            <root>
                <entry id="1" name="gene1" type="gene"/>
                <entry id="2" name="gene2" type="gene"/>
                <relation entry1="1" entry2="2" type="PPrel">
                    <subtype name="activation" value="--&gt;"/>
                    <subtype name="indirect effect" value="..&gt;"/>
                </relation>
            </root>
    '''
        xml_data = io.StringIO(xml_string)
        graph = KeggGraph(input_data=xml_data, type='gene', unique=False,
                                verbose=False)
        expected = pd.DataFrame([
            {'entry1': 'gene1', 'entry2': 'gene2', 'type': 'PPrel', 'subtype_name': 'activation', 'subtype_value': '-->',  'entry1_type': 'gene', 'entry2_type': 'gene'},
            {'entry1': 'gene1', 'entry2': 'gene2', 'type': 'PPrel', 'subtype_name': 'indirect effect', 'subtype_value': '..>', 'entry1_type': 'gene', 'entry2_type': 'gene'}])
        # assert pandas dataframe
        pd.testing.assert_frame_equal(graph.interaction, pd.DataFrame(expected))



    def test_parse_with_direction_ECrel(self):
        '''expect multiple rows for each relation if there are multiple subtypes'''
        # compound1 <-> gene1 <->  compound2 <-> gene2 -> compound3
        xml_string = '''
            <root>
                <entry id="1" name="compound1" type="compound"/>
                <entry id="2" name="gene1" type="gene"/>
                <entry id="3" name="compound2" type="compound"/>
                <entry id="4" name="gene2" type="gene"/>
                <entry id="5" name="compound3" type="compound"/>
                
                <relation entry1="2" entry2="4" type="ECrel">
                    <subtype name="compound" value="3"/>
                </relation>
                <reaction id="2" name="reaction2" type="reversible">
                    <substrate id="1"/>
                    <product id="3"/>
                </reaction>
                <reaction id="4" name="reaction4" type="reversible">
                    <substrate id="3"/>
                    <product id="5"/> 
                </reaction>
            </root>
    '''
        xml_data = io.StringIO(xml_string)
        graph = KeggGraph(input_data=xml_data, type='gene', unique=False,
                          verbose=False)
        expected = pd.DataFrame([
            {'entry1': 'gene1', 'entry2': 'gene2', 'type': 'PPrel', 'subtype_name': 'compound-propagation', 'subtype_value': 'custom', 'entry1_type': 'gene', 'entry2_type': 'gene'},
            {'entry1': 'gene2', 'entry2': 'gene1', 'type': 'PPrel', 'subtype_name': 'compound-propagation', 'subtype_value': 'custom', 'entry1_type': 'gene', 'entry2_type': 'gene'}])
        # assert pandas dataframe
        pd.testing.assert_frame_equal(graph.interaction, pd.DataFrame(expected).reset_index(drop=True))

    def test_parse_gene_only(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        output_df = parse_kgml_file(test_file, type="gene",
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_mpi(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        output_df = parse_kgml_file(test_file, type='mpi',
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_mpi_parser_directory(self):
        output_dir = Path('output')
        output_dir.mkdir(parents=True, exist_ok=True)
        runner = CliRunner()
        runner.invoke(network,
                      ['--type', 'mpi', 'hsa', test_file, '--out_dir', 'output', '--id_conversion', 'uniprot'])
        output_files = list(output_dir.glob('*.tsv'))
        assert len(output_files) > 0
        for file in output_files:
            df = pd.read_csv(file, sep='\t')
            assert not df.empty
            assert 'entry1' in df.columns
            assert 'entry2' in df.columns

    def test_get_edges_with_valid_data(self):
        parser = KeggGraph(input_data=self.xml_data, type='gene',
                           unique=False,
                           verbose=False)
        edges = parser.get_edges()
        assert not edges.empty
        assert 'entry1' in edges.columns
        assert 'entry2' in edges.columns
        assert 'type' in edges.columns
        assert 'subtype_name' in edges.columns
        assert 'subtype_value' in edges.columns

    def test_get_edges_with_empty_data(self):
        xml_data = '''
        <root></root>
        '''
        xml_data = io.StringIO(xml_data)
        with pytest.raises(FileNotFoundError):
            interaction = KeggGraph(input_data=xml_data, type='gene', unique=False,
                                    verbose=False)

    def test_get_edges_with_no_relations(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1" type="gene"/>
            <entry id="2" name="gene2" type="gene"/>
        </root>
        '''
        xml_data = io.StringIO(xml_data)
        with pytest.raises(FileNotFoundError):
            parser = KeggGraph(input_data=xml_data, type='gene', unique=False,
                               verbose=False)

    def test_remove_undefined_entries(self):
        output_df = parse_kgml_file(test_file, type="gene",
                                    unique=False, verbose=False)
        # check entries with undefined values
        assert output_df[output_df['entry1'].str.contains('undefined')].empty
        assert output_df[output_df['entry2'].str.contains('undefined')].empty

    def test_get_edges_with_mixed_data(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1" type="gene"/>
            <entry id="2" name="gene2" type="gene"/>
            <entry id="3" name="gene3" type="gene"/>
            <entry id="4" name="gene4" type="gene"/>
            <entry id="5" name="cpd:C00001" type="compound"/>
            <relation entry1="1" entry2="2" type="PPrel">
                <subtype name="compound" value="5"/>
            </relation>
            <relation entry1="3" entry2="4" type="PPrel">
                <subtype name="activation" value="act"/>
            </relation>
        </root>
        '''
        xml_data = io.StringIO(xml_data)
        parser = KeggGraph(input_data=xml_data, type='gene', unique=False,
                           verbose=False)
        edges = parser.get_edges()
        assert not edges.empty
        assert len(edges) == 2
        assert 'entry1' in edges.columns
        assert 'entry2' in edges.columns
        assert 'type' in edges.columns
        assert 'subtype_name' in edges.columns
        assert 'subtype_value' in edges.columns
        assert 'entry1_type' in edges.columns
        assert 'entry2_type' in edges.columns

    # test auto_relation_fix
    def test_auto_relation_fix(self):
        xml_data = '''
            <pathway name="path:hsa00010" org="hsa" number="00010"
                 title="Glycolysis / Gluconeogenesis"
                 image="https://www.kegg.jp/kegg/pathway/hsa/hsa00010.png"
                 link="https://www.kegg.jp/kegg-bin/show_pathway?hsa00010">
            <entry id="1" name="gene1" type="gene"/>
            <entry id="5" name="cpd:C00001" type="compound"/>
            <relation entry1="1" entry2="5" type="PPrel">
              <subtype name="activation" value="act"/>
            </relation>
        </pathway>
        '''
        xml_data = io.StringIO(xml_data)
        interaction = KeggGraph(input_data=xml_data, type='mpi', auto_correction='remove', unique=False,
                                id_converter=None, verbose=False)
        assert interaction.interaction.empty
        # interaction = Interaction(input_data=xml_data, type='mpi', auto_relation_fix='fix', unique=False, id_conversion=None, names=False, verbose=False)
        # assert len(interaction.data) == 1


# class TestGetMolEmbedding():
#
#     def test_get_smiles_gene_only(self):
#         parser = KeggGraph(input_data=test_file, type='gene')
#         smiles = parser.get_smiles()
#         assert smiles is None
#
#     def test_get_smiles_mpi(self):
#         parser = KeggGraph(input_data=test_file, type='mpi')
#         smiles = parser.get_smiles()
#         assert isinstance(smiles, dict)
#         assert len(smiles) > 0

    # def test_get_mol_embedding(self):
    #     graph = KeggGraph(input_data=test_file, type='mpi')
    #     graph.get_mol_embedding(transformer='morgan')
    #     assert isinstance(graph.mol_embedding, dict)
    #     assert len(graph.mol_embedding) > 0


# class TestProteinEmbedding:
#
#     def test_protein_embedding(self):
#         graph = KeggGraph(input_data=test_file, type='gene')
#         graph.get_protein_embedding(embedding_file='data/protein_embedding.h5')
#         assert isinstance(graph.protein_embedding, dict)
#         assert len(graph.protein_embedding) > 0


# class TestGetProteinFeatures:
#     def test_get_protein_features(self):
#         graph = KeggGraph(input_data=test_file, type='gene', id_converter=None)
#         # raise value error
#         with pytest.raises(ValueError):
#             graph.get_protein_features()

    # def test_get_protein_features_valid(self):
    #     converter = Converter('hsa', target='uniprot')
    #     graph = KeggGraph(input_data=test_file, type='gene', id_converter=converter)
    #     # expect  request.exception.HTTPError
    #     # with pytest.raises(requests.exceptions.HTTPError, match="500 Server Error"):
    #     graph.get_protein_features()
    #     assert isinstance(graph.protein_features, dict)
    #     assert len(graph.protein_features) > 0

    def test_gene_propagation(self):
        fake_instance = Mock()
        data = [
            ['gene1', 'compound1', 'PCrel', 'none', 'none', 'gene', 'compound'],
            ['compound1', 'compound2', 'PCrel', 'none', 'none', 'compound', 'compound'],
            ['compound2', 'gene3', 'PCrel', 'compound', 'none', 'compound', 'gene'],
        ]
        df = pd.DataFrame(data, columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        fake_instance.interaction = df
        result = KeggGraph.propagate(fake_instance, df, type_keep='gene')
        expected = pd.DataFrame([['gene1', 'gene3', 'PPrel', 'compound-propagation', 'custom', 'gene', 'gene']],
                                columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        # assert pandas dataframe
        pd.testing.assert_frame_equal(result, expected)


    def test_compound_propagation(self):
        fake_instance = Mock()
        data = [
            ['compound1', 'gene1', 'PCrel', 'none', 'none',  'compound', 'gene'],
            ['gene1', 'gene2', 'PPrel', 'none', 'none', 'gene', 'gene'],
            ['gene2', 'compound2', 'PCrel', 'compound', 'none', 'gene', 'compound'],
        ]
        df = pd.DataFrame(data, columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        fake_instance.interaction = df
        result = KeggGraph.propagate(fake_instance, df, type_keep='compound')
        expected = pd.DataFrame([['compound1', 'compound2', 'CCrel', 'gene-propagation', 'custom', 'compound', 'compound']],
                                columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        # assert pandas dataframe
        pd.testing.assert_frame_equal(result, expected)


def test_get_reactions():
    '''
    Test get_reactions method
    '''
    result = KeggGraph(input_data=test_file, type='gene').get_reactions()
    expected = pd.DataFrame([
        {'entry1': '176', 'entry2': '69', 'type': 'PCrel', 'name': 'reaction', 'value': 'rn:R04960', 'entry1_type': 'compound', 'entry2_type': 'gene'},
        {'entry1': '165', 'entry2': '69', 'type': 'PCrel', 'name': 'reaction', 'value': 'rn:R04960','entry1_type': 'compound', 'entry2_type': 'gene'},
        {'entry1': '69', 'entry2': '180', 'type': 'PCrel', 'name': 'reaction', 'value': 'rn:R04960', 'entry1_type': 'gene', 'entry2_type': 'compound'},
        {'entry1': '176', 'entry2': '165', 'type': 'CCrel', 'name': 'multi-substrate', 'value': 'rn:R04960', 'entry1_type': 'compound', 'entry2_type': 'compound'},
        {'entry1': '165', 'entry2': '176', 'type': 'CCrel', 'name': 'multi-substrate', 'value': 'rn:R04960', 'entry1_type': 'compound', 'entry2_type': 'compound'}])
    # assert pandas dataframe
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)


# def test_to_nextworkx():
#     result = KeggGraph(input_data=test_file, type='gene').to_networkx()
#     # plot the graph
#     assert result.number_of_nodes() == 24
