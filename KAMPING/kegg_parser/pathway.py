#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""

import re
import json
from typing import Union, Literal

import typer
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import urllib.request as request
from itertools import combinations
import xml.etree.ElementTree as ET
from collections import defaultdict
from KAMPING.kegg_parser import utils
from KAMPING.kegg_parser import convert
from KAMPING.kegg_parser import protein_metabolite_parser


class InteractionParser():
    """
   A class to parse gene interactions from KGML files.

   Attributes:
   ----------
   input_data : str
       Path to the input KGML file.
   mixed : bool
       Whether to include mixed interactions.
   unique : bool
       Whether to ensure unique interactions.
   id_conversion : Union[None, Literal['uniprot', 'ncbi']]
       Type of ID conversion to apply.
   names : bool
       Whether to include human-readable names.
   verbose : bool
       Whether to print verbose output.
   """

    def __init__(self,
                 input_data: str,
                 type: Literal['gene-only', 'MPI', 'original'],
                 unique: bool = False,
                 graphics: bool = False,
                 id_conversion: Union[None | Literal['uniprot', 'ncbi']] = None,
                 names: bool = False,
                 verbose: bool = False):
        """
      Initializes the GenesInteractionParser with the given parameters.

      Parameters:
      ----------
      input_data : str
          Path to the input KGML file
       type : Literal['gene-only', 'MPI', 'original']
            Type of interaction to parse
      unique : bool, optional
          Whether to ensure unique interactions (default is False).
      graphics : bool, optional
          Whether to include graphical information (default is False).
      id_conversion : Union[None, Literal['uniprot', 'ncbi']], optional
          Type of ID conversion to apply (default is None).
      names : bool, optional
          Whether to include human-readable names (default is False).
      verbose : bool, optional
          Whether to print verbose output (default is False).
      """

        self.id_conversion = id_conversion
        self.input_data = input_data
        self.type = type
        self.unique = unique
        self.names = names
        self.verbose = verbose

        tree = ET.parse(input_data)
        self.root = tree.getroot()

        self.conversion_dictionary = self._get_conversion_dictionary()
        if self.names:
            self.names_dictionary = self._get_names_dictionary(self.conversion_dictionary)


    def _get_edges(self):
        """
        Parses the KGML file to extract edges.

        Returns:
        -------
        pd.DataFrame
            DataFrame containing the edges.
        """
        pathway_link = self.root.get('link')

        d = []

        # Parse the relation and subtype elements
        # this will expand entry with subentry into multiple rows
        # for example
        # <relation entry1="20" entry2="37" type="PPrel">
        # <subtype name="activation" value="--&gt;"/>
        # <subtype name="indirect effect" value="..&gt;"/>
        # will be expanded into two rows
        # 20 37 PPrel activation --&gt;
        # 20 37 PPrel indirect effect ..&gt;
        for relation in self.root.findall('relation'):
            for subtype in relation:
                d1=relation.attrib
                d2=subtype.attrib
                d3=json.dumps(d1),json.dumps(d2)
                d.append(d3)

        edgelist=[]
        for line in d:
            x=line[0].replace("{","").replace("}","").replace('"entry1":',"") \
                .replace('"entry2":',"").replace('"type":',"").replace('"name":',"") \
                .replace('"value":',"").replace('"','').replace(',',"\t").replace(' ', '')
            y=line[1].replace("{","").replace("}","").replace('"name":',"") \
                .replace('"',"").replace('value:',"").replace(",",'\t').replace(' ', '')
            edgelist.append(x+"\t"+y)

        df=pd.DataFrame(edgelist)
        if df.empty:
            # throw error if no edges are found
            raise FileNotFoundError(f'ERROR: File "{self.input_data}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')

        df=df[0].str.split("\t", expand=True).rename({0: 'entry1',1: 'entry2',
                                                      2: 'types', 3:'name',
                                                      4: 'value'}, axis='columns')
        # reorder columns as entry1, entry2, types, value, name
        df = df[['entry1', 'entry2', 'types', 'value', 'name']]

        # convert compound value to kegg id if only relation.type is "compound"
        def apply_conversion(row):
            if row['name'] == 'compound':
                return self.conversion_dictionary.get(row['value'], row['value'])
            else:
                return row['value']
        df['value'] = df.apply(apply_conversion, axis=1)

        # Convert entry1 and entry2 id to kegg id
        df['entry1'] = df['entry1'].map(self.conversion_dictionary)
        df['entry2'] = df['entry2'].map(self.conversion_dictionary)
        # Split the entry1 and entry2 into lists
        # entry1 and entry2 can be a list of genes
        df['entry1'] = df['entry1'].astype(str).str.split(' ', expand = False)
        df['entry2'] = df['entry2'].astype(str).str.split(' ', expand = False)
        return df

    def _get_conversion_dictionary(self):
        if self.unique:
            conversion_dictionary = utils.conv_dict_unique(self.root)
        else:
            conversion_dictionary = utils.conv_dict(self.root)
        return conversion_dictionary

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


    def _add_names(self, df):
        '''
        This function adds the human-understandable names to the dataframe
        '''
        df['entry1_name'] = df.entry1.map(self.names_dictionary)
        df['entry2_name'] = df.entry2.map(self.names_dictionary)
        df.insert(1, 'entry1_name', df.pop('entry1_name'))
        df.insert(3, 'entry2_name', df.pop('entry2_name'))
        return df

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
        df = self._get_edges()

        # expode the entry1 and entry2 columns
        df_out = df.explode('entry1', ignore_index = True).explode('entry2', ignore_index = True)

        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not df_out[(df_out['entry1'].str.startswith('cpd:')) | (df_out['entry2'].str.startswith('cpd:')) | (df_out['entry1'].str.startswith('undefined')) | (df_out['entry2'].str.startswith('undefined'))].empty

        # if not mixed, remove "path" entries and propagate compounds
        if self.type == 'gene-only' or self.type == 'MPI':
            # Remove edges with "path" entries
            df_out = df_out[(~df_out['entry1'].str.startswith('path')) & (~df_out['entry2'].str.startswith('path'))]
            if self.type == "gene-only":
                if has_compounds_or_undefined:
                    df_out = self._propagate_compounds(df_out)
            if self.type == 'MPI':
                MPI_parser = protein_metabolite_parser.ProteinMetabliteParser(keep_PPI=True)
                df_out = MPI_parser.parse_dataframe(df_out)

        if self.id_conversion is not None:
            # convert the edges to the desired id type
            id_converter = convert.Converter(species=self.root.get('org'), target=self.id_conversion,
                                             unique=self.unique)
            df_out = id_converter.convert_dataframe(df_out)

        return df_out

def parse_pathway(input_data: str, wd: str,
                  type: Literal['gene-only', 'MPI', 'original'], unique: bool = False,
                  id_conversion: Union[None | Literal['unprot', 'ncbi']] = None,
                  names: bool = False, verbose: bool = False):
    '''
    Converts a folder of KGML files or a single KGML file into a weighted
    edgelist of genes that can be used in graph analysis.
    '''

    wd = Path(wd)
    # create the output directory if it does not exist
    wd.mkdir(parents=True, exist_ok=True)

    if Path(input_data).is_dir():
        for file in Path(input_data).glob('*.xml'):
            try:
                gip = InteractionParser(type=type, input_data=file,
                                        id_conversion=id_conversion,
                                        unique=unique, names=names,
                                        verbose=verbose)
                df_out = gip.parse_file()
                df_out.to_csv(wd / f'{file.stem}.tsv', sep='\t', index=False)
            except FileNotFoundError as e:
                typer.echo(typer.style(e, fg=typer.colors.RED, bold=True))
                continue
    else:
        gip = InteractionParser(type=type,
                                input_data=input_data,
                                id_conversion=id_conversion,
                                unique=unique, names=names,
                                verbose=verbose)
        df_out = gip.parse_file()
        df_out.to_csv(wd / f'{Path(input_data).stem}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    # for test
    output_df = parse_pathway('data/kgml_hsa/hsa00010.xml',
                              id_conversion='uniprot',
                              wd=Path.cwd(), mixed=True, unique=False, names=False, verbose=False)
    print(output_df)