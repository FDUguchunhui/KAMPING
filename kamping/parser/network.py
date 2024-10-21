#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""
import contextlib
import itertools
import logging
import os
from typing import Union, Any
import h5py
import requests.exceptions
from overrides import overrides
from typing_extensions import Literal
import typer
import pandas as pd
import networkx as nx
import xml.etree.ElementTree as ET
import torch_geometric as pyg
from rdkit import RDLogger
from kamping.uniprot.uniprot import submit_id_mapping, check_id_mapping_results_ready, get_id_mapping_results_link, \
    get_id_mapping_results_search

# Disable specific log
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')

from kamping.mol.utils import fetch_mol_file_string, get_smiles
from kamping.parser.utils import get_unique_compound_values, get_unique_proteins
from kamping.parser import utils
from kamping.parser import protein_metabolite_parser


class KeggGraph():
    ''''
    undefined nodes are removed from the final output
    '''

    def __init__(self,
                 input_data: str,
                 type: Literal['gene', 'metabolite', 'mpi'],
                 unique: bool = False,
                 protein_group_as_interaction: bool = True,
                 multi_substrate_as_interaction: bool = True,
                 auto_correction: Union[None, Literal['fix', 'remove']] = 'fix',
                 directed: bool = True,
                 id_converter: Any = None,
                 verbose: bool = False):
        '''
        Initialize the GenesInteractionParser object

        '''
        # todo: add function to filter common compounds
        self.auto_correction = auto_correction
        self.id_converter = id_converter
        self.input_data = input_data
        self.type = type
        self.unique = unique
        self.protein_group_as_interaction = protein_group_as_interaction
        self.multi_substrate_as_interaction = multi_substrate_as_interaction
        self.directed = directed
        self.verbose = verbose


        if self.verbose:
            typer.echo(typer.style(f"Now parsing: {self.root.get('title')}...", fg=typer.colors.GREEN, bold=False))
        self.root = ET.parse(input_data).getroot()
        self.entry_dictionary = utils.entry_id_conv_dict(self.root, unique=unique)
        self.edges = self.get_edges()
        self.mol_smiles = None
        self.mol_embedding = None
        self.protein_embedding = None

        if self.id_converter is not None:
            # convert the edges to the desired id type
            self.edges = self.id_converter.convert_dataframe(self.edges)

        # remove row with undefined entries, this maynot be necessary after group-expansion
        # self.edges = self.edges[~self.edges['entry1'].str.startswith('undefined')]
        # self.edges = self.edges[~self.edges['entry2'].str.startswith('undefined')]

        # auto_relation_fix
        if self.auto_correction is not None:
            self.auto_fix_relation()

        # remove prefix from entry1 and entry2
        self.edges['entry1'] = self.edges['entry1'].str.replace(r'^[^:]+:', '', regex=True)
        self.edges['entry2'] = self.edges['entry2'].str.replace(r'^[^:]+:', '', regex=True)

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
        # get type
        df['entry1_type'] = df['entry1'].map(lambda x: self.entry_dictionary.at[x, 'type'])
        df['entry2_type'] = df['entry2'].map(lambda x: self.entry_dictionary.at[x, 'type'])

        if df.empty:
            # throw error if no edges are found
            raise FileNotFoundError(f'ERROR: File "{self.input_data}" cannot be parsed.\nVisit {pathway_link} for pathway details.\nThere are likely no edges in which to parse...')
        # process the edge direction
        df = self.process_edge_direction(df)

        # remove interaction relationship with type "maplink"
        # replace ECrel with reaction step1: remove ECrel
        df_ECrel = df[df['type'] == 'ECrel']
        df = df[~df['type'].isin(['maplink', 'ECrel'])]

        # expand PCrel with reaction
        MPI_parser = protein_metabolite_parser.ProteinMetabliteParser(keep_PPI=True)
        df = MPI_parser.parse_dataframe(df)
        # replace ECrel with reaction step2: add reaction
        df = pd.concat([df, self.get_reactions()])
        if self.type == 'metabolite':
            # propagate the compound to get gene only interaction
            df = self.propagate(df, type_keep='compound')
        elif self.type == 'gene':
            df = self.propagate(df, type_keep='gene')

        # group expansion
        df = self.group_expansion(df)

        # convert entry id to kegg id
        df = self.convert_entry_id_to_kegg_id(df)

        # rename columns name to subtype_name and value to subtype_value
        df = df.rename(columns={'name': 'subtype_name', 'value': 'subtype_value'})

        # sort based on entry1 and entry2
        return df.sort_values(by=['entry1', 'entry2']).reset_index(drop=True)

    
    @property
    def interaction(self):
        # df = self.edges
        # # Identify rows where entry1 and entry2 values are swapped
        # df_sorted = df.apply(lambda row: tuple(sorted([row['entry1'], row['entry2'],row['type'], row['name'], row['value']])), axis=1)
        # df['sorted_entries'] = df_sorted
        # # Remove duplicate rows based on sorted entries
        # df = df.drop_duplicates(subset=['sorted_entries'])
        # # Set direction as "undirected" for rows with swapped values, otherwise "directed"
        # df.loc[:, 'direction'] = df_sorted.duplicated(keep=False).map({True: 'undirected', False: 'directed'})
        # return df.drop(columns=['sorted_entries']).reset_index(drop=True)
        return self.edges

    @property
    def proteins(self):
        # create node features
        return utils.get_unique_proteins(self.edges)

    @property
    def compounds(self):
        return  get_unique_compound_values(self.edges)

    @property
    def num_nodes(self):
        return len(self.proteins) + len(self.compounds)

    @property
    def num_edges(self):
        return len(self.edges)

    def __str__(self):
        # return the title of the graph and the number of nodes and edges
        return (f'''KEGG Pathway: 
            [Title]: {self.root.get('title')} \n 
            [Name]: {self.root.get('name')} \n 
            [Org]: {self.root.get('org')} \n 
            [Link]: {self.root.get('link')} \n 
            [Image]: {self.root.get('image')} \n 
            [Link]: {self.root.get('link')} \n 
            Number of Nodes: {self.num_nodes} \n 
            Number of Edges: {self.num_edges}''')


    def convert_entry_id_to_kegg_id(self, df: pd.DataFrame) -> pd.DataFrame:
        # convert compound value to kegg id if only relation.type is "compound"
        # compound is a list with one element mapped from dict
        df['value'] =df.apply(lambda row: self.entry_dictionary.at[row['value'], 'name'] if row['name'] == 'compound' else row['value'], axis=1)
        df['entry1'] = df['entry1'].map(lambda x: self.entry_dictionary.at[x, 'name'])
        df['entry2'] = df['entry2'].map(lambda x: self.entry_dictionary.at[x, 'name'])
        # expand entry1 and entry2 columns
        df = df.explode('entry1').explode('entry2').explode('value')
        return df

    def _get_names_dictionary(self, conversion_dictionary):
        '''
        Get the names dictionary for the given GenesInteractionParser object
        Returns a dictionary with the entry id as the key and the entry human-understandable name as the value.
        '''
        names_dictionary = utils.names_dict(self.root, self.root.get('org'), conversion_dictionary)
        return self.names_dictionary

    def propagate(self,df: pd.DataFrame,  type_keep: str) -> pd.DataFrame:
        '''
        Propagate the compound to get gene only interaction
        '''
        type_remove = None
        if type_keep == 'gene':
            type_remove = 'compound'
        elif type_keep == 'compound':
            type_remove = 'gene'

        G = nx.from_pandas_edgelist(df, source='entry1', target='entry2', edge_attr='name', create_using=nx.DiGraph())
        # create a dict with key in entry and value in entry_type
        entry1_type_dict = dict(zip(df['entry1'], df['entry1_type']))        # get all nodes with
        entry2_type_dict = dict(zip(df['entry2'], df['entry2_type']))
        entry_type_dict = {**entry1_type_dict, **entry2_type_dict}
        # set node attribute
        nx.set_node_attributes(G, entry_type_dict, 'entry_type')
        # get all nodes with type
        nodes_of_type = [node for node, data in G.nodes(data=True) if data['entry_type'] == type_keep]
        # use BFS to find all paths from gene to gene through GENE_PROPAGATION_TYPES
        new_edges = []
        for source in nodes_of_type:
            for target in nodes_of_type:
                if source != target:
                    # Find all paths from source to target
                    for path in nx.all_simple_paths(G, source=source, target=target):
                        # Check if all intermediate nodes are compounds
                        path_propagated_node_check = [G.nodes[node]['entry_type'] == type_remove for node in path[1:-1]]
                        # check if it not empty and all intermediate nodes are compounds
                        if path_propagated_node_check and all(path_propagated_node_check):
                            if type_keep == 'gene':
                                new_edges.append((source, target, 'PPrel', 'compound-propagation', 'custom', 'gene', 'gene'))
                            elif type_keep == 'compound':
                                new_edges.append((source, target, 'CCrel', 'gene-propagation', 'custom', 'compound', 'compound'))
        # add new edges to the df
        new_edges_df = pd.DataFrame(new_edges, columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        # drop duplicates based on entry1, entry2, type, name, value
        new_edges_df = new_edges_df.drop_duplicates(subset=['entry1', 'entry2', 'type', 'name', 'value'])
        df = pd.concat([df, new_edges_df], ignore_index=True).reset_index(drop=True)
        # reorder columns
        df = df[['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type']]
        # remove rows with entry1 or entry2 of type since they are original edges before propagation
        relation_with_type_remove = (df['entry1_type'] == type_remove) | (df['entry2_type'] == type_remove)
        df = df[~relation_with_type_remove].reset_index(drop=True)
        return df
    

    def auto_fix_relation(self):
        if self.type == 'gene-only':
            # when entry1_type == "compound" and entry2_type == "compound" remove the row
            self.edges = self.edges[~((self.edges['entry1_type'] == 'compound') & (self.edges['entry2_type'] == 'compound'))]
        if self.auto_correction == 'fix':
            # when entry1_type == "gene" and entry2_type == "gene" and type != "PPrel" or "GErel" then type = "PPrel"
            self.edges.loc[(self.edges['entry1_type'] == 'gene') & (self.edges['entry2_type'] == 'gene') & (self.edges['type'].isin(['PPrel', 'GErel', 'ECrel'])), 'type'] = 'PPrel'
            # when entry1_type == "gene" and entry2_type == "compound" and type != "PCrel" then type = "PCrel"
            self.edges.loc[(self.edges['entry1_type'] == 'gene') & (self.edges['entry2_type'] == 'compound') & (self.edges['type'] != 'PCrel'), 'type'] = 'PCrel'
            # when entry1_type == "compound" and entry2_type == "gene" and type != "PCrel" then type = "PCrel"
            self.edges.loc[(self.edges['entry1_type'] == 'compound') & (self.edges['entry2_type'] == 'gene') & (self.edges['type'] != 'PCrel'), 'type'] = 'PCrel'

        elif self.auto_correction == 'remove':
            # when entry1_type == "gene" and entry2_type == "gene" and type != "PPrel" or "GErel" then type = "PPrel"
            self.edges = self.edges[~(self.edges['entry1_type'] == 'gene') & (self.edges['entry2_type'] == 'gene') & (~self.edges['type'].isin(['PPrel', 'GErel', 'ECrel']))]
            # when entry1_type == "gene" and entry2_type == "compound" and type != "PCrel" then type = "PCrel"
            self.edges = self.edges[~(self.edges['entry1_type'] == 'gene') & (self.edges['entry2_type'] == 'compound') & (self.edges['type'] != 'PCrel')]
            # when entry1_type == "compound" and entry2_type == "gene" and type != "PCrel" then type = "PCrel"
            self.edges = self.edges[~(self.edges['entry1_type'] == 'compound') & (self.edges['entry2_type'] == 'gene') & (self.edges['type'] != 'PCrel')]



    def group_expansion(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        When there is entry with type "group", expand the group to the individual entries
        '''

        group_id_dict = utils.get_group_to_id_mapping(self.root)
        # for each value in group_id_dict, whici is a list of group members
        # create clique edges between the group members

        clique_edges = pd.DataFrame(columns=['entry1', 'entry2', 'type', 'name', 'value'])
        for group, members in group_id_dict.items():
            for member1, member2 in itertools.combinations(members, 2):
                # clique parsed as undirected represented by two directed edges
                clique_edges = pd.concat([clique_edges,
                                         pd.DataFrame.from_records([
                                             {'entry1': member1, 'entry2': member2, 'type': 'PPrel', 'name': 'protein-group', 'value': 'custom', 'entry1_type': 'gene', 'entry2_type': 'gene'},
                                             {'entry1': member2, 'entry2': member1, 'type': 'PPrel', 'name': 'protein-group', 'value': 'custom', 'entry1_type': 'gene', 'entry2_type': 'gene'}])])
        df['entry1'] =   [group_id_dict[entry]  if entry in group_id_dict else entry for entry in df['entry1']]
        df['entry2'] =   [group_id_dict[entry]  if entry in group_id_dict else entry for entry in df['entry2']]
        df = df.explode('entry1').explode('entry2')
        # toggle group-clique
        if self.protein_group_as_interaction:
            df = pd.concat([df, clique_edges])
        # remove row with entry1_type="group" and entry2_type="group"
        df = df[~((df['entry1_type'] == 'group') | (df['entry2_type'] == 'group'))]
        return df

    def process_edge_direction(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Process the direction of the edges
        '''
        # for type=PPrel or GErel, name="state change" or "binding/association", "dissociation", "missing interaction",
        # set the direction to undirected otherwise set the direction to directed
        undirected_PPrel = ['compound', 'state change', 'binding/association', 'dissociation', 'missing interaction']
        df_PPrel_reverse = df.loc[(df['type'].isin(['PPrel', 'GErel', 'PCrel'])) & (df['name'].isin(undirected_PPrel))]
        # reverse the entry1 and entry2 value
        df_PPrel_reverse.assign(entry1=df_PPrel_reverse['entry2'], entry2=df_PPrel_reverse['entry1'])
        # add the reversed edge to the original df
        df = pd.concat([df, df_PPrel_reverse])
        return df


    def get_reactions(self):
        '''
        Get the reaction from the graph

          <reaction id="426" name="rn:R04355" type="irreversible">
        <substrate id="425" name="cpd:C03939"/>
        <substrate id="366" name="cpd:C01209"/>
        <product id="367" name="cpd:C05744"/>
        '''
        reactions = reactions = [
            {
                'id': reaction.attrib['id'],
                'name': reaction.attrib['name'],
                'type': reaction.attrib['type'],
                'substrate': [substrate.attrib['id'] for substrate in reaction.findall('substrate')],
                'product': [product.attrib['id'] for product in reaction.findall('product')]
            }
            for reaction in self.root.findall('reaction')
        ]

        reactions = pd.DataFrame(reactions, columns=['id', 'name', 'type', 'substrate', 'product'])
        # add column direction
        reactions['direction'] = reactions['type'].map(lambda x: 'undirected' if x == 'reversible' else 'directed')
        # expand each row in reactions into two rows substrate -> id, and id -> product
        new_edges = []
        for index, row in reactions.iterrows():
            new_edges.append({'entry1': row['substrate'], 'entry2': row['id'], 'type': 'PCrel',
                              'name': 'reaction', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'gene'})
            new_edges.append({'entry1': row['id'], 'entry2': row['product'], 'type': 'PCrel',
                              'name': 'reaction', 'value': row['name'], 'entry1_type': 'gene', 'entry2_type': 'compound'})
            if row['type'] == 'reversible':
                new_edges.append({'entry1': row['id'], 'entry2': row['substrate'], 'type': 'PCrel',
                                  'name': 'reaction', 'value': row['name'], 'entry1_type': 'gene', 'entry2_type': 'compound'})
                new_edges.append({'entry1': row['product'], 'entry2': row['id'], 'type': 'PCrel',
                              'name': 'reaction', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'gene'})
            # add cliques edges between substrate
            if self.multi_substrate_as_interaction:
                if len(row['substrate']) > 1:
                    for substrate1, substrate2 in itertools.combinations(row['substrate'], 2):
                        new_edges.append({'entry1': substrate1, 'entry2': substrate2, 'type': 'CCrel',
                                          'name': 'multi-substrate', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'compound'})
                        new_edges.append({'entry1': substrate2, 'entry2': substrate1, 'type': 'CCrel',
                                          'name': 'multi-substrate', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'compound'})

                if len(row['product']) > 1:
                    for product1, product2 in itertools.combinations(row['product'], 2):
                        new_edges.append({'entry1': product1, 'entry2': product2, 'type': 'CCrel',
                                          'name': 'multi-product', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'compound'})
                        new_edges.append({'entry1': product2, 'entry2': product1, 'type': 'CCrel',
                                          'name': 'multi-product', 'value': row['name'], 'entry1_type': 'compound', 'entry2_type': 'compound'})
        new_edges = pd.DataFrame(new_edges, columns=['entry1', 'entry2', 'type', 'name', 'value',  'entry1_type', 'entry2_type'])
        # explode entry1 and entry2
        new_edges = new_edges.explode('entry1').explode('entry2')
        # convert entry1 and entry2 to kegg id
        # if direction is undirected expand it into two rows with reversed entry1 and entry2 and direction as undirected


        return new_edges

    def save_mol_embedding(self, file_path):
        '''
        Save the molecular embeddings to a file
        '''
        if self.mol_embedding is None:
            raise ValueError('Molecular embeddings have not been generated yet!')
        # create parent directory if not exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        # save  mol_ids, and mol_embedding as h5
        with h5py.File(file_path, 'w') as h5file:
            for key, value in self.embeddings.items():
                h5file.create_dataset(key, data=value)




    def get_protein_features(self) -> dict[str, str]:
        '''
        Get the protein features
        '''
        if self.id_converter is None or self.id_converter.target != 'uniprot':
            raise ValueError('The target must be set to uniprot')

        proteins = get_unique_proteins(self.edges)
        # when 500 error the server is down
        try:
            job_id = submit_id_mapping(
                from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=proteins
            )
        except requests.exceptions.HTTPError as e:
            typer.echo(typer.style(f'Error when submitting job: {e}. '
                                   f'Please https://www.uniprot.org/id-mapping see if server is down.', fg=typer.colors.RED, bold=True))
            raise e

        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)

        queries = []
        sequences = []

        # logging failed ids
        if 'failedIds' in results:
            logging.warning(results['failedIds'])
        for result in results['results']:
            queries.append(result['from'])
            # code below return None if the key is not found
            sequence = result.get('to', {}).get('sequence', {}).get('value', '')
            sequences.append(sequence)

        self.protein_features = dict(zip(queries, sequences))


    def to_networkx(self) -> nx.DiGraph:
        '''
        Convert the graph to networkx
        '''
        G = nx.from_pandas_edgelist(self.expand_undirected_edges(),
                                    source='entry1', target='entry2', edge_attr=['name', 'type', 'value', 'direction', 'entry1_type', 'entry2_type',], create_using=nx.DiGraph())

        # set node attribute
        # create dict with key in self.proteins and value as "gene"

        # create dict with key in self.molecules and value as "compound"
        # combine the two dictionaries
        node_attributes = {**{protein: 'gene' for protein in self.proteins},
                           **{molecule: 'compound' for molecule in self.compounds}}

        nx.set_node_attributes(G, node_attributes, name='type')
        return G


    def expand_undirected_edges(self) -> pd.DataFrame:
        '''
        Expand the unidirectional edges
        '''
        # if indirection is "undirected", add the reverse edge
        return pd.concat([self.edges, self.edges[self.edges['direction'] == 'undirected'].rename(columns={'entry1': 'entry2', 'entry2': 'entry1'})])


    def to_homogenous_pyg_data(self):
        '''
        Convert the graph to PyG data
        '''
        if self.type == 'gene':
            # convert the graph to networkx
            G = self.to_networkx()
            # add embedding to the node attributes
            nx.set_node_attributes(G, self.protein_embedding, name='embedding')
            # convert the graph to PyG data
            data =  pyg.utils.from_networkx(G, node_attrs=['embedding'])
        elif self.type == 'metabolite':
            # convert the graph to networkx
            G = self.to_networkx()
            # add embedding to the node attributes
            nx.set_node_attributes(G, self.mol_embedding, name='embedding')
            # convert the graph to PyG data
            data =  pyg.utils.from_networkx(G, node_attrs=['embedding'])

        return data




    #
    # def load_edge_csv(self, data: pd.DataFrame, protein_mapping:dict, mol_mapping:str):
    #
    #     # for df_pc if entry2 is a protein, swap entry1 and entry2
    #     # todo: find a more elegant way to do this
    #     # add it to parser function
    #     # this will remove direction information while it is impossible to get direction from the KGML file
    #     data.loc[(data['type'] == 'PCrel') & (data['entry2_type'] == 'gene'), ['entry1', 'entry2']] = df.loc[(df['type'] == 'PCrel') & (df['entry2'].str.startswith('up')), ['entry2', 'entry1']].values
    #
    #
    # # Remove rows with entry1 or entry2 not in the mapping
    # mapping_keys = list(protein_mapping.keys()) + list(mol_mapping.keys())
    #
    # # remove prefix from entry1 and entry2
    # df['entry1'] = df['entry1'].str.replace(r'^[^:]+:', '', regex=True)
    # df['entry2'] = df['entry2'].str.replace(r'^[^:]+:', '', regex=True)
    #
    # # remove rows with entry1 or entry2 not in the mapping
    # df = df[df['entry1'].isin(mapping_keys) & df['entry2'].isin(mapping_keys)]
    #
    #
    # # split df into two parts: one for protein-protein edges and the other for protein-compound edges
    # # the PP edges will have type "PPrel"
    # # the PC edges will have type "PCrel"
    # df_pp = df[df['type'] == 'PPrel']
    # df_pc = df[df['type'] == 'PCrel']
    #
    #
    # # df_pc edges
    # pc_src = [protein_mapping[entry] for entry in df_pc['entry1']]
    # pc_dst = [mol_mapping[entry] for entry in df_pc['entry2']]
    # pc_index = torch.tensor([pc_src, pc_dst])
    #
    # # df_pp edges
    # pp_src = [protein_mapping[entry] for entry in df_pp['entry1']]
    # pp_dst = [protein_mapping[entry] for entry in df_pp['entry2']]
    # pp_index = torch.tensor([pp_src, pp_dst])
    #
    # pc_edge_attr = None
    # pp_edge_attr = None
    #
    # return pc_index, pc_edge_attr, pp_index, pp_edge_attr