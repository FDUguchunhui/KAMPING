#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""

import json
from typing import Union
from typing_extensions import Literal
import typer
import pandas as pd
import networkx as nx
import xml.etree.ElementTree as ET
from kamping.parser import utils
from kamping.parser import convert
from kamping.parser import protein_metabolite_parser


class InteractionParser():
    ''''
    undefined nodes are removed from the final output
    '''

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
        df['value'] =df.apply(lambda row: self.conversion_dictionary.at[row['value'], 'name'] if row['name'] == 'compound' else row['value'], axis=1)

        # add another two columns to indicate type of entry1 and entry2
        df['entry1_type'] = df['entry1'].map(lambda x: self.conversion_dictionary.at[x, 'type'])
        df['entry2_type'] = df['entry2'].map(lambda x: self.conversion_dictionary.at[x, 'type'])
        df['entry1'] = df['entry1'].map(lambda x: self.conversion_dictionary.at[x, 'name'])
        df['entry2'] = df['entry2'].map(lambda x: self.conversion_dictionary.at[x, 'name'])



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

        # perform last clean-up to make list of  [cpd: ...] into a single string
        # those are mapped by the conversion dictionary for column entry1, entry2, and value
        df = df.explode('entry1').explode('entry2').explode('value')


        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not df[(df['entry1'].str.startswith('cpd:')) | (df['entry2'].str.startswith('cpd:')) | (df['entry1'].str.startswith('undefined')) | (df['entry2'].str.startswith('undefined'))].empty

        # if not mixed, remove "path" entries and propagate compounds
        if self.type == 'gene-only' :
            # Remove edges with "path" entries
            df = df[(~df['entry1'].str.startswith('path')) & (~df['entry2'].str.startswith('path'))]
            if has_compounds_or_undefined:
                df = self._propagate_compounds(df)
        elif self.type == 'MPI':
            # remove interaction relationship with type "maplink"
            # which will leave "ECrel", "GErel", and "PPrel", and "PCrel"
            df = df[df['type'] != 'maplink']
            MPI_parser = protein_metabolite_parser.ProteinMetabliteParser(keep_PPI=True)
            df = MPI_parser.parse_dataframe(df)
        elif self.type == 'original':
            pass
        else:
            raise ValueError(f'Invalid type: {self.type}')

        if self.id_conversion is not None:
            # convert the edges to the desired id type
            id_converter = convert.Converter(species=self.root.get('org'), target=self.id_conversion,
                                             unique=self.unique)
            df = id_converter.convert_dataframe(df)

        # remove row with undefined entries
        df = df[~df['entry1'].str.startswith('undefined')]
        df = df[~df['entry2'].str.startswith('undefined')]

        return df


