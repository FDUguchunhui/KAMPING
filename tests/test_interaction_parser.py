import pytest
import pandas as pd
from pathlib import Path
from KAMPING.kegg_parser.pathway import InteractionParser, parse_pathway


def parse_kgml_file(file_path, **kwargs):
    parser = InteractionParser(input_data=file_path, **kwargs)
    return parser.parse_file()


class TestGenesInteractionParser:

    test_file = 'data/hsa00010_test.xml'

    def test_parse_empty_kgml_file(self):
        with pytest.raises(FileNotFoundError):
            parse_kgml_file('data/kgml_hsa/empty.xml', type='gene-only')

    def test_genes_parser_single_file(self):
        output_df = parse_kgml_file(self.test_file, type="gene-only", unique=False, verbose=False)
        assert not output_df.empty
        assert 'entry1' in output_df.columns
        assert 'entry2' in output_df.columns

    def test_genes_parser_directory(self):
        output_dir = Path('output')
        output_dir.mkdir(parents=True, exist_ok=True)
        parse_pathway(self.test_file, wd='output', type="gene-only", unique=False, verbose=False)
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
        parse_pathway(self.test_file, wd='output', type="MPI", unique=False,  verbose=False)
        output_files = list(output_dir.glob('*.tsv'))
        assert len(output_files) > 0
        for file in output_files:
            df = pd.read_csv(file, sep='\t')
            assert not df.empty
            assert 'entry1' in df.columns
            assert 'entry2' in df.columns