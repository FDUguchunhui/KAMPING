#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for converting pathway TSV files into UniProt and NCBI IDs
"""

import json
import pathlib
import re
from pathlib import Path

import pandas as pd
import typer
from pandas import DataFrame
from typing import Literal, overload
from KAMPING.kegg_parser.utils import get_conversion_dictionary



pd.options.mode.chained_assignment = None
app = typer.Typer()




class Converter:
    def __init__(self, species: str,
                 target: Literal['uniprot', 'ncbi'] = 'uniprot',
                 unique: bool = False,
                 unmatched: Literal['drop', 'keep'] = 'drop',
                 verbose: bool = False):
        # todo: folder or file should be input for function
        self.species = species
        self.unique = unique
        self.target = target
        self.unmatched = unmatched
        self.conversion = get_conversion_dictionary(self.species, target=target)

    def _process_dataframe(self, df):
        if self.unique:
            # Extract the terminal modifiers and create a new column
            # This enables the re-addition of the modifiers at the
            # end of the function.
            df['match1'] = df['entry1'].str.extract(r'(-[0-9]+)')
            df['match2'] = df['entry2'].str.extract(r'(-[0-9]+)')
            # Remove the terminal modifier so that the IDs map properly
            # to the KEGG API call
            df['entry1'] = df['entry1'].str.replace(r'(-[0-9]+)', '', regex=True)
            df['entry2'] = df['entry2'].str.replace(r'(-[0-9]+)', '', regex=True)

        # Map to convert KEGG IDs to target IDs. Note lists are returned
        # for some conversions.
        df['entry1_conv'] = df['entry1'].map(self.conversion)
        df['entry2_conv'] = df['entry2'].map(self.conversion)

        # umatched machanism
        if self.unmatched == 'drop':
            # drop rows with unmatched gene entries
            df = df[~((df['entry1'].str.startswith('hsa')) & (df['entry1_conv'].isna()))]
            df = df[~((df['entry2'].str.startswith('hsa')) & (df['entry2_conv'].isna()))]
            df['entry1'] = df['entry1_conv'].fillna(df['entry1'])
            df['entry2'] = df['entry2_conv'].fillna(df['entry2'])
        elif self.unmatched == 'keep':
            # Fills nans with entries from original columns
            df['entry1'] = df['entry1_conv'].fillna(df['entry1'])
            df['entry2'] = df['entry2_conv'].fillna(df['entry2'])


        # Drop the extra column as it's all now in entry1/2 columns
        df = df.drop(['entry1_conv', 'entry2_conv'], axis=1)

        # Due to one to many mapping, we need to explode the lists
        df = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)
        df['entry1'] = df['entry1'].str.replace(f'{self.target}:', '')
        df['entry2'] = df['entry2'].str.replace(f'{self.target}:', '')

        if self.unique:
            df['entry1'] = df['entry1'] + df['match1']
            df['entry2'] = df['entry2'] + df['match2']
            df = df.drop(['match1', 'match2'], axis=1)

        return df

    def convert_dataframe(self, df: DataFrame) -> DataFrame:
        '''
        Converts a dataframe of KEGG IDs to UniProt or NCBI IDs
        '''
        df_out = self._process_dataframe(df)
        return df_out

    def convert_file(self, input_data: str)  -> DataFrame:
        '''
        Wrapper function for converting a single file
        '''
        file = Path(input_data)
        df = pd.read_csv(input_data, delimiter='\t')
        typer.echo(f"Now converting {file.name} to {self.target} IDs...")
        df_out = self._process_dataframe(df)
        return df_out

        # df_out.to_csv(out_dir / 'up_{}'.format(file.name), sep='\t', index=False)
        # typer.echo(f'Now converting {file.name} to NCBI IDs...')
            # df_out.to_csv(self.wd / 'ncbi_{}'.format(file.name), sep='\t', index=False)

        # print work done
        typer.echo(typer.style(f'Conversion of {file.name} complete!', fg=typer.colors.GREEN, bold=True))


def genes_convert(species: str,
                  input_data: str,
                  wd: str,
                  target: Literal['uniprot', 'ncbi'] = 'uniprot',
                  unique: bool = False,
                  verbose: bool = False):
    '''
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis.
    '''
    # set working directory
    wd = Path(wd)
    wd.mkdir(parents=True, exist_ok=True)

    converter = Converter(species, unique=unique,
                          target=target,
                          verbose=verbose)

    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.tsv'):
            try:
                df_converted = converter.convert_file(file)
                df_converted.to_csv(Path(wd) / f'{file.stem}_{target}.tsv', sep='\t', index=False)
            except FileNotFoundError as e:
                typer.echo(typer.style(e, fg=typer.colors.RED, bold=True))
                continue
    else:
        df_converted = converter.convert_file(input_data)


if __name__ == '__main__':
    # set working directory
    # for test
    genes_convert('hsa',
                  input_data='data/mixed_output_PC/hsa00010_mpi.tsv',
                  wd='data/mixed_output_converted',
                  target='uniprot')