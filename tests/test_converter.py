import pandas as pd
import pytest
from kamping.parser.convert import Converter
from kamping.parser.network import KeggGraph


class TestConverter:
    def test_convert_uniprot(self, snapshot):
        converter = Converter(species='hsa', target='uniprot')
        graph = KeggGraph(input_data='data/hsa_test.xml', type='gene',
        unique=False,
        verbose=False)
        converter.convert(graph)
        assert graph.interaction is not None
        snapshot.assert_match(graph.interaction)

    def test_convert_ncbi(self, snapshot):
        converter = Converter(species='hsa', target='ncbi')
        graph = KeggGraph(input_data='data/hsa_test.xml', type='gene',
        unique=False,
        verbose=False)
        converter.convert(graph)
        # assert graph.interaction is not None
        assert graph.interaction is not None
        snapshot.assert_match(graph.interaction)

