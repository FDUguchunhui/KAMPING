import os
import pandas as pd
import pytest
import pandas.testing as pdt

from kamping.parser.utils import entry_id_conv_dict, get_group_to_id_mapping
from kamping.utils import read_all_tsv_files
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
        expected = pd.DataFrame.from_dict( {'1': {'name': ['gene1-1', 'gene2-1'], 'type':'gene'}, '2': {'name':['gene3-2'], 'type':'gene'}}, orient='index')
        # set index name
        expected.index.name = 'id'
        pdt.assert_frame_equal(result, expected)


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
        expected = pd.DataFrame.from_dict( {'1': {'name': ['gene1', 'gene2'], 'type':'gene'}, '2': {'name':['gene3'], 'type':'gene'}}, orient='index')
       # set index name
        expected.index.name = 'id'
        pdt.assert_frame_equal(result, expected)

    def test_entry_id_conv_dict_with_empty_entries(self):
        xml_data = '''
        <root></root>
        '''
        root = ET.fromstring(xml_data)
        result = entry_id_conv_dict(root, unique=True)
        assert result.empty

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

        expected = pd.DataFrame.from_dict( {'1': {'name': ['gene1-1', 'gene2-1'], 'type':'gene'}, '2': {'name':['gene3-2'], 'type':'gene'},
                                            '3': {'name':['compound1-3'], 'type':'compound'}}, orient='index')
    # set index name
        expected.index.name = 'id'
        pdt.assert_frame_equal(result, expected)


def test_get_group_to_id_mapping():
    xml = '''
    <pathway>
            <entry id="352" name="undefined" type="group">
        <graphics fgcolor="#000000" bgcolor="#FFFFFF"
             type="rectangle" x="804" y="1206" width="92" height="51"/>
        <component id="240"/>
        <component id="241"/>
        <component id="242"/>
        <component id="243"/>
        <component id="244"/>
    </entry>
    </pathway>
    '''
    root = ET.fromstring(xml)
    result = get_group_to_id_mapping(root)
    expected = {"352": ["240", "241", "242", "243", "244"]}
    assert result == expected


