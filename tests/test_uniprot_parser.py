import pytest
from unittest.mock import patch
from kamping.uniprot.parser import get_protein_seq, get_protein_seqs

class TestUniprotParser:
    def test_get_protein_seq_with_valid_accession(self):
        sequence = get_protein_seq("Q96C23")
        assert sequence == ('MASVTRAVFGELPSGGGTVEKFQLQSDLLRVDIISWGCTITALEVKDRQGRASDVVLGFA'
                            'ELEGYLQKQPYFGAVIGRVANRIAKGTFKVDGKEYHLAINKEPNSLHGGVRGFDKVLWTP'
                            'RVLSNGVQFSRISPDGEEGYPGELKVWVTYTLDGGELIVNYRAQASQATPVNLTNHSYFN'
                            'LAGQASPNINDHEVTIEADTYLPVDETLIPTGEVAPVQGTAFDLRKPVELGKHLQDFHLN'
                            'GFDHNFCLKGSKEKHFCARVHHAASGRVLEVYTTQPGVQFYTGNFLDGTLKGKNGAVYPK'
                            'HSGFCLETQNWPDAVNQPRFPPVLLRPGEEYDHTTWFKFSVA')
