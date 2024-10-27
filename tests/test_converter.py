import pandas as pd
import pytest
from kamping.parser.convert import Converter
from kamping.parser.network import KeggGraph
from kamping.parser.convert import conv_dict_from_url


class TestConverter:
    def test_convert_uniprot(self, snapshot):
        converter = Converter(species='hsa', gene_target='uniprot')
        graph = KeggGraph(input_data='data/hsa_test.xml', type='gene',
        unique=False,
        verbose=False)
        converter.convert(graph)
        assert graph.interaction is not None
        snapshot.assert_match(graph.interaction)

    def test_convert_ncbi(self, snapshot):
        converter = Converter(species='hsa', gene_target='ncbi', compound_target='pubchem')
        graph = KeggGraph(input_data='data/hsa_test.xml', type='mixed',
        unique=False,
        verbose=False)
        converter.convert(graph)
        # assert graph.interaction is not None
        assert graph.interaction is not None
        snapshot.assert_match(graph.interaction)


    def test_from_response_to_dict_gene(self):
        url = 'http://rest.kegg.jp/conv/hsa/ncbi-geneid'
        conv_dict = conv_dict_from_url(url)
        assert conv_dict is not None


    def test_from_response_to_dict_compound(self):
        url = 'http://rest.kegg.jp/conv/compound/pubchem'
        conv_dict = conv_dict_from_url(url)
        assert conv_dict is not None

    def test_get_compound_conv_dict(self):
        converter = Converter(species='hsa', gene_target='uniprot', compound_target='pubchem')
        compound_conv_dict = converter.get_compound_conv_dict()
        assert compound_conv_dict is not None