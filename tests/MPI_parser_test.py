import unittest
import pandas as pd
from src.parser.MPI_parser import MPIParser, expand_EC_relation


class MPIParserTest(unittest.TestCase):
    def test_expand_EC_relation(self):
        # create pd.Dataframe
        # entry1 entry2 type value name
        #  hsa:130589	hsa:2538	ECrel	cpd:C00267-90	compound
        # create pd.Dataframe from the above data
        input = pd.DataFrame.from_dict({'entry1': ['hsa:130589'],
                                        'entry2': ['hsa:2538'],
                      'type': ['ECrel'], 'value': ['cpd:C00267-90'], 'name': ['compound']})
        # write test case for expand_EC_relation
        expected = pd.DataFrame.from_dict({
            'entry1': ['hsa:130589', 'hsa:2538'],
             'entry2': ['cpd:C00267-90', 'cpd:C00267-90'],
             'type': ['ECrel', 'ECrel'],
             'value': ['cpd:C00267-90', 'cpd:C00267-90'],
             'name': ['compound', 'compound']}
        )

        output = expand_EC_relation(input)
        # assert output and expected are equal and they are pd.Dataframe
        self.assertTrue(expected.equals(output))

if __name__ == '__main__':
    unittest.main()
