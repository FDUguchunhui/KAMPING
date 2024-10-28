#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for converting pathway TSV files into UniProt and NCBI IDs
"""

import json
import logging
import pathlib
import re
from collections import defaultdict
from pathlib import Path
from urllib import request as request

import pandas as pd
import typer
from pandas import DataFrame
from typing import Literal, overload, Union

pd.options.mode.chained_assignment = None
app = typer.Typer()



class Converter:
    ''''
    Class to convert KEGG IDs to other IDs

    Parameters
    ----------
    species: str
        The species to convert the KEGG IDs to
    gene_target: str
        The target gene ID to convert to. Options are 'uniprot', 'ncbi', 'kegg', default is 'uniprot'
    compound_target: str
        The target compound ID to convert to. Options are 'pubchem', 'chebi', 'kegg', default is 'kegg'
    unique: bool
        Whether to keep the terminal modifiers in the IDs
    unmatched: str
        What to do with unmatched gene entries. Options are 'drop' or 'keep'
    verbose: bool
        Whether to print out the progress of the conversion

    Usage
    -----
    converter = Converter(species='hsa', gene_target='uniprot', compound_target='pubchem')
    converter.convert(graph)

    '''
    def __init__(self, species: str,
                 gene_target: Literal['uniprot', 'ncbi', 'kegg'] = 'kegg',
                 compound_target: Literal['pubchem', 'chebi', 'kegg'] = 'kegg',
                 unique: bool = False,
                 unmatched: Literal['drop', 'keep'] = 'drop',
                 verbose: bool = True):
        self.logger = logging.getLogger(__name__)
        if verbose:
            self.logger.setLevel(logging.INFO)
        self.species = species
        self.unique = unique
        self.gene_target = gene_target
        self.compound_target = compound_target
        self.unmatched = unmatched

        self.gene_conv_dict = None
        if self.gene_target != 'kegg':
            self.gene_conv_dict = self.get_gene_conv_dict()
        self.compound_conv_dict = None
        # convert compouund ID if asked
        if self.compound_target != 'kegg':
            # append compound conversion dict to gene conversion dict
            self.compound_conv_dict = self.get_compound_conv_dict()

    def _process_dataframe(self, graph):
        df = graph.edges
        graph_name = graph.root.get('name')

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
        if self.gene_conv_dict is not None:
            df.loc[df['entry1_type'] == 'gene', 'entry1_conv'] = df.loc[df['entry1_type'] == 'gene', 'entry1'].map(self.gene_conv_dict)
            df.loc[df['entry2_type'] == 'gene', 'entry2_conv'] = df.loc[df['entry2_type'] == 'gene', 'entry2'].map(self.gene_conv_dict)
        else:
            # if no gene conversion is needed, fill with original values
            df.loc[df['entry1_type'] == 'gene', 'entry1_conv'] = df.loc[df['entry1_type'] == 'gene', 'entry1']
            df.loc[df['entry2_type'] == 'gene', 'entry2_conv'] = df.loc[df['entry2_type'] == 'gene', 'entry2']

        if self.compound_conv_dict is not None:
            df.loc[df['entry1_type'] == 'compound', 'entry1_conv'] = df.loc[df['entry1_type'] == 'compound', 'entry1'].map(self.compound_conv_dict)
            df.loc[df['entry2_type'] == 'compound', 'entry2_conv'] = df.loc[df['entry2_type'] == 'compound', 'entry2'].map(self.compound_conv_dict)
        else:
            # if no compound conversion is needed, fill with original values
            df.loc[df['entry1_type'] == 'compound', 'entry1_conv'] = df.loc[df['entry1_type'] == 'compound', 'entry1']
            df.loc[df['entry2_type'] == 'compound', 'entry2_conv'] = df.loc[df['entry2_type'] == 'compound', 'entry2']



        # umatched machanism
        if self.unmatched == 'drop':
            # drop rows with unmatched gene entries
            # NAN conversion are present in form of empty list
            df = df[~(df['entry1_conv'].isna())]
            df = df[~(df['entry2_conv'].isna())]
            df['entry1'] = df['entry1_conv']
            df['entry2'] = df['entry2_conv']

        elif self.unmatched == 'keep':
                # Fills nans with entries from original columns
                df['entry1'] = df['entry1_conv'].fillna(df['entry1'])
                df['entry2'] = df['entry2_conv'].fillna(df['entry2'])

        # Drop the extra column as it's all now in entry1/2 columns
        df = df.drop(['entry1_conv', 'entry2_conv'], axis=1)

        # Due to one to many mapping, we need to explode the lists
        df = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)

        if self.unique:
            df['entry1'] = df['entry1'] + df['match1']
            df['entry2'] = df['entry2'] + df['match2']
            df = df.drop(['match1', 'match2'], axis=1)

        # check if there are any NA entries
        if (df[['entry1', 'entry2']].isna().any().any()):
            # remove row with entry1 or entry NA and print warning
            df = df.dropna(subset=['entry1', 'entry2'])
            self.logger.warning(f'{graph_name} dropped due to NA entries')

        # log work done
        self.logger.info(f'Conversion of {graph_name} complete!')
        return df


    def convert_file(self, input_data: str)  -> DataFrame:
        '''
        Wrapper function for converting a single file
        '''
        file = Path(input_data)
        df = pd.read_csv(input_data, delimiter='\t')
        typer.echo(f"Now converting {file.name} to {self.gene_target} IDs...")
        df_out = self._process_dataframe(df)
        return df_out

        # df_out.to_csv(out_dir / 'up_{}'.format(file.name), sep='\t', index=False)
        # typer.echo(f'Now converting {file.name} to NCBI IDs...')
            # df_out.to_csv(self.wd / 'ncbi_{}'.format(file.name), sep='\t', index=False)

        # print work done
        typer.echo(typer.style(f'Conversion of {file.name} complete!', fg=typer.colors.GREEN, bold=True))

    def convert(self, graph) -> None:
        '''
        Converts a graph of KEGG IDs to UniProt or NCBI IDs
        '''
        if graph.gene_id_type == self.gene_target and graph.compound_id_type == self.compound_target:
            logging.info(f'''{graph.root.get("name")} already in gene: {self.gene_target} format
                         and compound: {self.compound_target} format. No nothing has been done.''')
            return
        graph.edges = self._process_dataframe(graph)
        graph.gene_id_type = self.gene_target
        graph.compound_id_type = self.compound_target


    def get_gene_conv_dict(self):
        '''
        Convert KEGG gene IDs to either NCBI gene IDs or UniProt IDs.
        '''
        if self.gene_target == 'uniprot':
            url = 'http://rest.kegg.jp/conv/%s/uniprot'
        elif self.gene_target == 'ncbi':
            url = 'http://rest.kegg.jp/conv/%s/ncbi-geneid'
        return conv_dict_from_url(url % self.species)


    def get_compound_conv_dict(self):
        '''
        Convert KEGG compound IDs to PubChem compound IDs.
        '''
        url_compound = f'http://rest.kegg.jp/conv/compound/{self.compound_target}'
        url_glycan = f'http://rest.kegg.jp/conv/glycan/{self.compound_target}'
        url_drug = f'http://rest.kegg.jp/conv/drug/{self.compound_target}'
        compound_conv_dict = {**conv_dict_from_url(url_compound), **conv_dict_from_url(url_glycan), **conv_dict_from_url(url_drug)}
        return compound_conv_dict


def conv_dict_from_url(url):
    response = request.urlopen(url).read().decode('utf-8')
    response = response.rstrip().rsplit('\n')
    kegg = []
    target = []
    for resp in response:
        target_id, kegg_id = resp.rsplit()
        target.append(target_id)
        kegg.append(kegg_id)

    d = dict()
    for key, value in zip(kegg, target):
        # initialize the key if it is not in the dictionary
        if key not in d:
            d[key] = []
        d[key].append(value)
    return d