import pytest
import pandas as pd
from pathlib import Path
from KAMPING.kegg_parser.protein_metabolite_parser import ProteinMetabliteParser, expand_relation_ECrel, expand_relation_PPrel, remove_suffix

def test_parse_dataframe_with_empty_input():
    parser = ProteinMetabliteParser()
    df = pd.DataFrame(columns=['entry1', 'entry2', 'types', 'value', 'name'])
    with pytest.raises(ValueError, match='input dataframe does not contain data to parse!'):
        parser.parse_dataframe(df)

def test_parse_dataframe_with_valid_input():
    parser = ProteinMetabliteParser()
    df = pd.DataFrame({
        'entry1': ['hsa:130589', 'hsa:2538'],
        'entry2': ['cpd:C00267-90', 'hsa:2538'],
        'types': ['ECrel', 'PPrel'],
        'value': ['cpd:C00267-90', 'cpd:C00267-90'],
        'name': ['compound', 'compound']
    })
    result = parser.parse_dataframe(df)
    assert not result.empty
    assert 'entry1' in result.columns
    assert 'entry2' in result.columns

def test_parse_file_with_nonexistent_file():
    parser = ProteinMetabliteParser()
    with pytest.raises(FileNotFoundError, match='File nonexistent_file.tsv not found!'):
        parser.parse_file('nonexistent_file.tsv', '.')

def test_parse_file_with_valid_file(tmp_path):
    parser = ProteinMetabliteParser()
    df = pd.DataFrame({
        'entry1': ['hsa:130589', 'hsa:2538'],
        'entry2': ['cpd:C00267-90', 'hsa:2538'],
        'types': ['ECrel', 'PPrel'],
        'value': ['cpd:C00267-90', 'cpd:C00267-90'],
        'name': ['compound', 'compound']
    })
    file_path = tmp_path / 'valid_file.tsv'
    df.to_csv(file_path, sep='\t', index=False)
    parser.parse_file(file_path, tmp_path)
    output_file = tmp_path / 'valid_file_mpi.tsv'
    assert output_file.exists()
    result = pd.read_csv(output_file, sep='\t')
    assert not result.empty
    assert 'entry1' in result.columns
    assert 'entry2' in result.columns

def test_expand_relation_ECrel_with_valid_input():
    input = pd.Series({
        'entry1': 'hsa:130589',
        'entry2': 'hsa:2538',
        'types': 'ECrel',
        'value': 'cpd:C00267-90',
        'name': 'compound'
    })
    expected = pd.DataFrame.from_dict({
        'entry1': ['hsa:130589', 'hsa:2538'],
        'entry2': ['cpd:C00267-90', 'cpd:C00267-90'],
        'types': ['PCrel', 'PCrel'],
        'value': ['custom', 'custom'],
        'name': ['enzyme-enzyme expansion', 'enzyme-enzyme expansion']
    })
    output = expand_relation_ECrel(input)
    pd.testing.assert_frame_equal(expected, output)

def test_expand_relation_PPrel_with_valid_input():
    input = pd.Series({
        'entry1': 'hsa:130589',
        'entry2': 'hsa:2538',
        'types': 'PPrel',
        'value': 'cpd:C00267-90',
        'name': 'compound'
    })
    expected = pd.DataFrame.from_dict({
        'entry1': ['hsa:130589', 'hsa:2538'],
        'entry2': ['cpd:C00267-90', 'cpd:C00267-90'],
        'types': ['PCrel', 'PCrel'],
        'value': ['custom', 'custom'],
        'name': ['protein-protein expansion', 'protein-protein expansion']
    })
    output = expand_relation_PPrel(input)
    pd.testing.assert_frame_equal(expected, output)

def test_remove_suffix_from_entry():
    assert remove_suffix('hsa:130589-90') == 'hsa:130589'
    assert remove_suffix('cpd:C00267-90') == 'cpd:C00267'
    assert remove_suffix('gl:12345-90') == 'gl:12345'