#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""

import json
from typing import Union, Literal

import typer
import pandas as pd
import networkx as nx
import lxml.etree as ET
from KAMPING.kegg_parser import utils
from KAMPING.kegg_parser import convert
from KAMPING.kegg_parser import protein_metabolite_parser


class InteractionParser():

    def __init__(self,
                 input_data: str,
                 type: Literal['gene-only', 'MPI', 'original'],
                 unique: bool = False,
                 id_conversion: Union[Literal['uniprot', 'ncbi'], None] = None,
                 names: bool = False,
                 verbose: bool = False):
        '''
        Initialize the GenesInteractionParser object

        '''

        self.id_conversion = id_conversion
        self.input_data = input_data
        self.type = type
        self.unique = unique
        self.names = names
        self.verbose = verbose

        tree = ET.parse(input_data)
        self.root = tree.getroot()

        self.conversion_dictionary = utils.entry_id_conv_dict(self.root, unique=unique)

    def get_edges(self) -> pd.DataFrame:
        """
        Parses the KGML file to extract edges.

        Returns:
        -------
        pd.DataFrame
            DataFrame containing the edges.
        """
        pathway_link = self.root.get('link')

        # Parse the relation and subtype elements
        relations = [
            {
                **relation.attrib,
                **subtype.attrib
            }
            for relation in self.root.findall('relation')
            for subtype in relation
        ]

        # Create DataFrame from parsed relations
        df = pd.DataFrame(relations, columns=['entry1', 'entry2', 'type', 'name', 'value'])

        if df.empty:
            # throw error if no edges are found
            raise FileNotFoundError(f'ERROR: File "{self.input_data}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')

        # convert compound value to kegg id if only relation.type is "compound"
        # compound is a list with one element mapped from dict
        df['value'] = df.apply(lambda row: self.conversion_dictionary.get(row['value']) if row['name'] == 'compound' else row['value'], axis=1)

        # Convert entry1 and entry2 id to kegg id
        df['entry1'] = df['entry1'].map(self.conversion_dictionary)
        df['entry2'] = df['entry2'].map(self.conversion_dictionary)

        # expand entry1 and entry2 columns
        df = df.explode('entry1').explode('entry2')

        return df


    def _get_names_dictionary(self, conversion_dictionary):
        '''
        Get the names dictionary for the given GenesInteractionParser object
        Returns a dictionary with the entry id as the key and the entry human-understandable name as the value.
        '''
        names_dictionary = utils.names_dict(self.root, self.root.get('org'), conversion_dictionary)
        return self.names_dictionary




    # def _parse_multi_to_multi(self, df):
    #     '''
    #     This function takes a dataframe with multiple genes in each entry1 and entry2
    #     and returns a dataframe with each gene pair as a separate row
    #     '''
    #     # todo: use df explode instead
    #     edges = []
    #     for index, row in df.iterrows():
    #         entry1_list = row.iloc[0]
    #         entry2_list = row.iloc[1]
    #         for entry1 in entry1_list:
    #             for entry2 in entry2_list:
    #                 edges.append([entry1, entry2]  + row[2:].tolist())
    #
    #     df_out = pd.DataFrame(edges, columns = ['entry1', 'entry2', 'type', 'value', 'name'])
    #     return df_out

    def _propagate_compounds(self, xdf):
        G = nx.from_pandas_edgelist(xdf, source='entry1', target='entry2', edge_attr='name', create_using=nx.DiGraph())
        new_edges = []
        for node in G.nodes:
            if node.startswith(('cpd', 'undefined')) and node not in [n for n, d in G.in_degree() if d == 0] + [n for n, d in G.out_degree() if d == 0]:
                for i in G.in_edges(node):
                    for o in G.out_edges(node):
                        if not any(x.startswith(('cpd', 'undefined', 'path')) for x in [i[0], o[1]]):
                            new_edges.append([i[0], o[1], 'CPp', 'Custom', 'compound propagation'])
                        else:
                            for root in [n for n, d in G.in_degree() if d == 0]:
                                for leaf in [n for n, d in G.out_degree() if d == 0]:
                                    if nx.has_path(G, root, node) and nx.has_path(G, node, leaf):
                                        rpath, lpath = nx.shortest_path(G, root, node), nx.shortest_path(G, node, leaf)
                                        if not all(x.startswith(('cpd', 'undefined', 'path')) for x in rpath + lpath):
                                            rindex, lindex = [i for i, x in enumerate(rpath) if not x.startswith(('cpd', 'undefined', 'path'))], [i for i, x in enumerate(lpath) if not x.startswith(('cpd', 'undefined', 'path'))]
                                            new_edges.append([rpath[max(rindex)], lpath[min(lindex)], 'CPp', 'Custom', 'compound propagation'])
        df0 = pd.concat([xdf, pd.DataFrame(new_edges, columns=['entry1', 'entry2', 'type', 'value', 'name'])]).drop_duplicates()
        return df0[~df0['entry1'].str.startswith(('cpd', 'undefined', 'path')) & ~df0['entry2'].str.startswith(('cpd', 'undefined', 'path'))]


    def parse_file(self) -> pd.DataFrame:
        '''
        This function parses the KGML file and returns a dataframe of the edges

        Parameters:
        input_data: str
            The path to the KGML file
        out_dir: str
            The path to the output directory
        unique: bool
            If True, the output will contain unique nodes

        Returns: pd.DataFrame
        '''
        title = self.root.get('title')
        pathway = self.root.get('name').replace('path:', '')
        pathway_link = self.root.get('link')

        # Common operations
        if self.verbose:
            typer.echo(typer.style(f'Now parsing: {title}...', fg=typer.colors.GREEN, bold=False))
        df = self.get_edges()

        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not df[(df['entry1'].str.startswith('cpd:')) | (df['entry2'].str.startswith('cpd:')) | (df['entry1'].str.startswith('undefined')) | (df['entry2'].str.startswith('undefined'))].empty

        # if not mixed, remove "path" entries and propagate compounds
        if self.type == 'gene-only' :
            # Remove edges with "path" entries
            df = df[(~df['entry1'].str.startswith('path')) & (~df['entry2'].str.startswith('path'))]
            if self.type == "gene-only":
                if has_compounds_or_undefined:
                    df = self._propagate_compounds(df)
        if self.type == 'MPI':
            MPI_parser = protein_metabolite_parser.ProteinMetabliteParser(keep_PPI=True)
            df = MPI_parser.parse_dataframe(df)
        elif self.type == 'original':
            pass
        else:
            raise ValueError(f'Invalid type: {self.type}')


        # perform last clean-up to make list of  [cpd: ...] into a single string
        df = df.explode('entry1').explode('entry2')

        if self.id_conversion is not None:
            # convert the edges to the desired id type
            id_converter = convert.Converter(species=self.root.get('org'), target=self.id_conversion,
                                             unique=self.unique)
            df = id_converter.convert_dataframe(df)

        #todo: remove Undefined nodes


        return df


