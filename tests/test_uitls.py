import os
import pandas as pd
import pytest

from KAMPING.kegg_parser.utils import entry_id_conv_dict
from KAMPING.utils import read_all_tsv_files
import xml.etree.ElementTree as ET

def test_read_all_tsv_files_empty_directory(tmp_path):
    os.makedirs(tmp_path / 'empty_dir')
    result = read_all_tsv_files(tmp_path / 'empty_dir')
    assert result.empty

def test_read_all_tsv_files_single_file(tmp_path):
    df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
    file_path = tmp_path / 'file1.tsv'
    df.to_csv(file_path, sep='\t', index=False)
    result = read_all_tsv_files(tmp_path)
    assert not result.empty
    assert 'source' in result.columns
    assert result['source'].iloc[0] == 'file1'

def test_read_all_tsv_files_multiple_files(tmp_path):
    df1 = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
    df2 = pd.DataFrame({'col1': [5, 6], 'col2': [7, 8]})
    file_path1 = tmp_path / 'file1.tsv'
    file_path2 = tmp_path / 'file2.tsv'
    df1.to_csv(file_path1, sep='\t', index=False)
    df2.to_csv(file_path2, sep='\t', index=False)
    result = read_all_tsv_files(tmp_path)
    assert not result.empty
    assert 'source' in result.columns

def read_all_tsv_files_ignore_hidden_files(tmp_path):
    df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
    file_path = tmp_path / '.hidden_file.tsv'
    df.to_csv(file_path, sep='\t', index=False)
    result = read_all_tsv_files(tmp_path)
    assert result.empty



class TestEntryIdConvDict:
    def test_entry_id_conv_dict_with_unique_entries(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1 gene2" type="gene"/>
            <entry id="2" name="gene3" type="gene"/>
        </root>
        '''
        root = ET.fromstring(xml_data)
        result = entry_id_conv_dict(root, unique=True)
        assert result == {'1': ['gene1-1', 'gene2-1'], '2': ['gene3-2']}


    def test_entry_id_conv_dict_with_non_unique_entries(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1 gene2" type="gene"/>
            <entry id="2" name="gene3" type="gene"/>
        </root>
        '''
        root = ET.fromstring(xml_data)
        # non-unique entries
        result = entry_id_conv_dict(root, unique=False)
        assert result == {'1': ['gene1', 'gene2'], '2': ['gene3']}

    def test_entry_id_conv_dict_with_empty_entries(self):
        xml_data = '''
        <root></root>
        '''
        root = ET.fromstring(xml_data)
        result = entry_id_conv_dict(root, unique=True)
        assert result == {}

    def test_entry_id_conv_dict_with_mixed_entries(self):
        xml_data = '''
        <root>
            <entry id="1" name="gene1 gene2" type="gene"/>
            <entry id="2" name="gene3" type="gene"/>
            <entry id="3" name="compound1" type="compound"/>
        </root>
        '''
        root = ET.fromstring(xml_data)
        result = entry_id_conv_dict(root, unique=True)
        assert result == {'1': ['gene1-1', 'gene2-1'], '2': ['gene3-2'], '3': ['compound1-3']}