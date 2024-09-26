import unittest
import pandas as pd
from src.parser.protein_metabolite_parser import expand_relation_ECrel, expand_relation_PPrel
from pandas.testing import assert_frame_equal

class MPIParserTest(unittest.TestCase):
    def test_expand_relation_ECrel(self):
        # create pd.Dataframe
        # entry1 entry2 type value name
        #  hsa:130589	hsa:2538	ECrel	cpd:C00267-90	compound
        # create pd.Dataframe from the above data
        input = pd.Series({'entry1': 'hsa:130589',
                                        'entry2': 'hsa:2538',
                      'type': 'ECrel', 'value': 'cpd:C00267-90', 'name': 'compound'})
        # write test case for expand_EC_relation

        expected = pd.DataFrame.from_dict({
            'entry1': ['hsa:130589', 'hsa:2538'],
            'entry2': ['cpd:C00267-90', 'cpd:C00267-90'],
            'type': ['PCrel', 'PCrel'],
            'value': ['custom', 'custom'],
            'name': ['enzyme-enzyme expansion', 'enzyme-enzyme expansion']}
        )

        output = expand_relation_ECrel(input)
        # assert output and expected are equal and they are pd.Dataframe
        assert_frame_equal(expected, output)

    # test for expand_relation_PPrel
    def test_expand_relation_PPrel(self):
        # create pd.Dataframe
        # entry1 entry2 type value name
        #  hsa:130589	hsa:2538	ECrel	cpd:C00267-90	compound
        # create pd.Dataframe from the above data
        input = pd.Series({'entry1': 'hsa:130589',
                                        'entry2': 'hsa:2538',
                      'type': 'PPrel', 'value': 'cpd:C00267-90', 'name': 'compound'})
        # write test case for expand_EC_relation

        expected = pd.DataFrame.from_dict({
            'entry1': ['hsa:130589', 'hsa:2538'],
            'entry2': ['cpd:C00267-90', 'cpd:C00267-90'],
            'type': ['PCrel', 'PCrel'],
            'value': ['custom', 'custom'],
            'name': ['protein-protein expansion', 'protein-protein expansion']}
        )

        output = expand_relation_PPrel(input)
        # assert output and expected are equal and they are pd.Dataframe
        assert_frame_equal(expected, output)

if __name__ == '__main__':
    unittest.main()
