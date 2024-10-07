import io

import pytest
import pandas as pd
from pathlib import Path
from kamping.parser.network import InteractionParser
from kamping.main import network
import xml.etree.ElementTree as ET
from kamping.main import  Type, Identifier_Conversion

def parse_kgml_file(file_path, **kwargs):
    parser = InteractionParser(input_data=file_path, **kwargs)
    return parser.parse_file()


class TestGenesInteractionParser:

    test_file = 'data/hsa00010_test.xml'

    xml_data = '''
    <root>
        <relation entry1="1" entry2="2">
            <subtype name="compound" value="cpd:C00001"/>
        </relation>
    </root>
    '''
    # xml_data to file-like object
    xml_data = io.StringIO(xml_data)

    def test_parse_empty_kgml_file(self):
        with pytest.raises(OSError):
            parse_kgml_file('data/kgml_hsa/empty.xml', type='gene-only')

    def test_genes_parser_single_file(self):
        output_df = parse_kgml_file(self.test_file, type="gene-only", unique=False, verbose=False)
        assert not output_df.empty
        assert 'entry1' in output_df.columns
        assert 'entry2' in output_df.columns

    def test_genes_parser_directory(self):
        output_dir = Path('output')
        output_dir.mkdir(parents=True, exist_ok=True)
        network('gene-only', self.test_file, out_dir='output', id_conversion=None, unique=False, verbose=False)
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
        output_df = parse_kgml_file(self.test_file, type="MPI", unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_with_id_conversion(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        output_df = parse_kgml_file(self.test_file, type="gene-only", id_conversion='uniprot',
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_gene_only(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        output_df = parse_kgml_file(self.test_file, type="gene-only", id_conversion='uniprot',
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_parse_MPI(self, snapshot):
        '''
        Test with id_conversion='uniprot'
        '''
        output_df = parse_kgml_file(self.test_file, type='MPI', id_conversion='uniprot',
                                    unique=False, verbose=False)
        snapshot.assert_match(output_df)

    def test_MPI_parser_directory(self):
        output_dir = Path('output')
        output_dir.mkdir(parents=True, exist_ok=True)
        network("MPI", self.test_file, out_dir='output', id_conversion='uniprot', unique=False, verbose=False)
        output_files = list(output_dir.glob('*.tsv'))
        assert len(output_files) > 0
        for file in output_files:
            df = pd.read_csv(file, sep='\t')
            assert not df.empty
            assert 'entry1' in df.columns
            assert 'entry2' in df.columns

    def test_get_edges_with_valid_data(self):


        parser = InteractionParser(input_data=self.xml_data, type='gene-only', unique=False, id_conversion=None, names=False, verbose=False)
        edges = parser.get_edges()
        assert not edges.empty
        assert 'entry1' in edges.columns
        assert 'entry2' in edges.columns
        assert 'type' in edges.columns
        assert 'name' in edges.columns
        assert 'value' in edges.columns

    def test_get_edges_with_empty_data(self):
        xml_data = '''
        <root></root>
        '''
        xml_data = io.StringIO(xml_data)
        parser = InteractionParser(input_data=xml_data, type='gene-only', unique=False, id_conversion=None, names=False, verbose=False)
        with pytest.raises(FileNotFoundError):
            parser.get_edges()

    def test_get_edges_with_no_relations(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1" type="gene"/>
            <entry id="2" name="gene2" type="gene"/>
        </root>
        '''
        xml_data = io.StringIO(xml_data)
        parser = InteractionParser(input_data=xml_data, type='gene-only', unique=False, id_conversion=None, names=False, verbose=False)
        with pytest.raises(FileNotFoundError):
            parser.get_edges()

    def test_remove_undefined_entries(self):
        output_df = parse_kgml_file(self.test_file, type="gene-only", id_conversion='uniprot',
                                    unique=False, verbose=False)
        # check entries with undefined values
        assert output_df[output_df['entry1'].str.contains('undefined')].empty
        assert output_df[output_df['entry2'].str.contains('undefined')].empty

    def test_get_edges_with_mixed_data(self):
        xml_data = '''
        <root>
            <relation entry1="1" entry2="2">
                <subtype name="compound" value="cpd:C00001"/>
            </relation>
            <relation entry1="3" entry2="4">
                <subtype name="activation" value="act"/>
            </relation>
        </root>
        '''
        xml_data = io.StringIO(xml_data)
        parser = InteractionParser(input_data=xml_data, type='gene-only', unique=False, id_conversion=None, names=False, verbose=False)
        edges = parser.get_edges()
        assert not edges.empty
        assert len(edges) == 2
        assert 'entry1' in edges.columns
        assert 'entry2' in edges.columns
        assert 'type' in edges.columns
        assert 'name' in edges.columns
        assert 'value' in edges.columns




