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
import numpy as np
import requests.exceptions
import scikit_mol.fingerprints
from rdkit.Chem import PandasTools
from typing_extensions import Literal
import typer
import pandas as pd
import networkx as nx
import xml.etree.ElementTree as ET

from rdkit import RDLogger

from kamping.parser.utils import load_embedding_from_h5, get_unique_proteins
from kamping.uniprot.uniprot import submit_id_mapping, check_id_mapping_results_ready, get_id_mapping_results_link, \
    get_id_mapping_results_search

# Disable specific log
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')

from kamping.mol.utils import get_unique_compound_values, fetch_mol_file_string, get_smiles
from kamping.parser import utils
from kamping.parser import protein_metabolite_parser

# KGML types
# for relation, only gene, compound, group, map
# for species-specfic 'enzyme', 'reaction', 'reaction', 'brite', 'other' is not included
ENTRY_TYPES = [ 'ortholog', 'gene', 'group', 'compound', 'map']
# ignore  'maplink' for now
RELATION_TYPES = ['PPrel', 'GErel', 'PCrel', 'ECrel']
GENE_PROPAGATION_TYPES = ['compound', 'map']


class KeggGraph():
    ''''
    undefined nodes are removed from the final output
    '''

    def __init__(self,
                 input_data: str,
                 type: Literal['gene-only', 'mpi', 'original'],
                 auto_relation_fix: Union[None, Literal['fix', 'remove']] = 'fix',
                 unique: bool = False,
                 id_converter: Any = None,
                 names: bool = False,
                 verbose: bool = False):
        '''
        Initialize the GenesInteractionParser object

        '''
        self.auto_relation_fix = auto_relation_fix
        self.id_converter = id_converter
        self.input_data = input_data
        self.type = type
        self.unique = unique
        self.names = names
        self.verbose = verbose


        if self.verbose:
            typer.echo(typer.style(f"Now parsing: {self.root.get('title')}...", fg=typer.colors.GREEN, bold=False))
        self.root = ET.parse(input_data).getroot()
        self.graph_title = self.root.get('title')
        self.entry_dictionary = utils.entry_id_conv_dict(self.root, unique=unique)
        self.interaction = self.get_edges()
        self.mol_smiles = None
        self.mol_embedding = None
        self.protein_embedding = None

        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not self.interaction[(self.interaction['entry1'].str.startswith('cpd:')) | (self.interaction['entry2'].str.startswith('cpd:')) | (self.interaction['entry1'].str.startswith('undefined')) | (self.interaction['entry2'].str.startswith('undefined'))].empty

        # if not mixed, remove "path" entries and propagate compounds
        if self.type == 'gene-only' :
            # Remove edges with "path" entries
            self.interaction = self.interaction[self.interaction['type'] != 'maplink']
            if has_compounds_or_undefined:
                self.interaction = self.propagate_to_gene(self.interaction)
        elif self.type == 'mpi':
            # remove interaction relationship with type "maplink"
            # which will leave "ECrel", "GErel", and "PPrel", and "PCrel"
            self.interaction = self.interaction[self.interaction['type'] != 'maplink']
            MPI_parser = protein_metabolite_parser.ProteinMetabliteParser(keep_PPI=True)
            self.interaction = MPI_parser.parse_dataframe(self.interaction)
        elif self.type == 'original':
            pass
        else:
            raise ValueError(f'Invalid type: {self.type}')

        if self.id_converter is not None:
            # convert the edges to the desired id type
            self.interaction = self.id_converter.convert_dataframe(self.interaction)

        # remove row with undefined entries, this maynot be necessary after group-expansion
        # self.interaction = self.interaction[~self.interaction['entry1'].str.startswith('undefined')]
        # self.interaction = self.interaction[~self.interaction['entry2'].str.startswith('undefined')]

        # auto_relation_fix
        if self.auto_relation_fix is not None:
            self.auto_fix_relation()

        # remove prefix from entry1 and entry2
        self.interaction['entry1'] = self.interaction['entry1'].str.replace(r'^[^:]+:', '', regex=True)
        self.interaction['entry2'] = self.interaction['entry2'].str.replace(r'^[^:]+:', '', regex=True)

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


        df = self.group_expansion(df)
        df = self.convert_entry_id_to_kegg_id(df)


        return df

    def convert_entry_id_to_kegg_id(self, df: pd.DataFrame) -> pd.DataFrame:
        # convert compound value to kegg id if only relation.type is "compound"
        # compound is a list with one element mapped from dict
        df['value'] =df.apply(lambda row: self.entry_dictionary.at[row['value'], 'name'] if row['name'] == 'compound' else row['value'], axis=1)
        # add another two columns to indicate type of entry1 and entry2
        df['entry1_type'] = df['entry1'].map(lambda x: self.entry_dictionary.at[x, 'type'])
        df['entry2_type'] = df['entry2'].map(lambda x: self.entry_dictionary.at[x, 'type'])
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

    def propagate_to_gene(self):
        '''
        Propagate the compound to get gene only interaction
        '''
        G = nx.from_pandas_edgelist(self.interaction, source='entry1', target='entry2', edge_attr='name', create_using=nx.DiGraph())
        # create a dict with key in entry and value in entry_type
        entry1_type_dict = dict(zip(self.interaction['entry1'], self.interaction['entry1_type']))        # get all nodes with
        entry2_type_dict = dict(zip(self.interaction['entry2'], self.interaction['entry2_type']))
        entry_type_dict = {**entry1_type_dict, **entry2_type_dict}
        # set node attribute
        nx.set_node_attributes(G, entry_type_dict, 'entry_type')
        # get all nodes with type "compound" or "path"
        # get all nodes with type "gene"
        genes = [node for node, data in G.nodes(data=True) if data['entry_type'] == 'gene']
        # use BFS to find all paths from gene to gene through GENE_PROPAGATION_TYPES
        new_edges = []
        for source in genes:
            for target in genes:
                if source != target:
                    # Find all paths from source to target
                    for path in nx.all_simple_paths(G, source=source, target=target):
                        # Check if all intermediate nodes are compounds
                        if all(G.nodes[node]['entry_type'] == 'compound' for node in path[1:-1]):
                            new_edges.append((source, target, 'PPrel', 'compound-propagation', 'custom', 'gene', 'gene'))
        # add new edges to the self.interaction
        new_edges_df = pd.DataFrame(new_edges, columns=['entry1', 'entry2', 'type', 'name', 'value', 'entry1_type', 'entry2_type'])
        self.interaction = pd.concat([self.interaction, new_edges_df], ignore_index=True).reset_index(drop=True)
        # remove rows with entry1 and entry2 with type "compound"
        relation_with_compound = (self.interaction['entry1_type'] == 'compound') | (self.interaction['entry2_type'] == 'compound')
        self.interaction = self.interaction[~relation_with_compound].reset_index(drop=True)


    def auto_fix_relation(self):
        if self.auto_relation_fix == 'fix':
            # when entry1_type == "gene" and entry2_type == "gene" and type != "PPrel" or "GErel" then type = "PPrel"
            self.interaction.loc[(self.interaction['entry1_type'] == 'gene') & (self.interaction['entry2_type'] == 'gene') & (~self.interaction['type'].isin(['PPrel', 'GErel', 'ECrel']))] = 'PPrel'
            # when entry1_type == "gene" and entry2_type == "compound" and type != "PCrel" then type = "PCrel"
            self.interaction.loc[(self.interaction['entry1_type'] == 'gene') & (self.interaction['entry2_type'] == 'compound') & (self.interaction['type'] != 'PCrel')] = 'PCrel'
            # when entry1_type == "compound" and entry2_type == "gene" and type != "PCrel" then type = "PCrel"
            self.interaction.loc[(self.interaction['entry1_type'] == 'compound') & (self.interaction['entry2_type'] == 'gene') & (self.interaction['type'] != 'PCrel')] = 'PCrel'
            # when entry1_type == "compound" and entry2_type == "compound" remove the row
            # todo: need to mute this when try to get metablolite only
            self.interaction = self.interaction[~((self.interaction['entry1_type'] == 'compound') & (self.interaction['entry2_type'] == 'compound'))]
        elif self.auto_relation_fix == 'remove':
            # when entry1_type == "gene" and entry2_type == "gene" and type != "PPrel" or "GErel" then type = "PPrel"
            self.interaction = self.interaction[~(self.interaction['entry1_type'] == 'gene') & (self.interaction['entry2_type'] == 'gene') & (~self.interaction['type'].isin(['PPrel', 'GErel', 'ECrel']))]
            # when entry1_type == "gene" and entry2_type == "compound" and type != "PCrel" then type = "PCrel"
            self.interaction = self.interaction[~(self.interaction['entry1_type'] == 'gene') & (self.interaction['entry2_type'] == 'compound') & (self.interaction['type'] != 'PCrel')]
            # when entry1_type == "compound" and entry2_type == "gene" and type != "PCrel" then type = "PCrel"
            self.interaction = self.interaction[~(self.interaction['entry1_type'] == 'compound') & (self.interaction['entry2_type'] == 'gene') & (self.interaction['type'] != 'PCrel')]
            # when entry1_type == "compound" and entry2_type == "compound" remove the row
            # todo: need to mute this when try to get metablolite only
            self.interaction = self.interaction[~((self.interaction['entry1_type'] == 'compound') & (self.interaction['entry2_type'] == 'compound'))]


    def get_mol_embedding(self, transformer, dim=1024, **kwargs) -> dict[list, np.array]:
        '''
        Get the molecular embeddings
        '''
        if self.mol_embedding is None:
            if transformer == 'morgan':
                transformer = scikit_mol.fingerprints.MorganFingerprintTransformer(nBits=dim, **kwargs)
            elif transformer == 'rdkit':
                transformer = scikit_mol.fingerprints.RDKitFingerprintTransformer(fpSize=dim, **kwargs)
            elif transformer == 'atom-pair':
                transformer = scikit_mol.fingerprints.AtomPairFingerprintTransformer(nBits=dim, **kwargs)
            elif transformer == 'topological':
                transformer = scikit_mol.fingerprints.TopologicalTorsionFingerprintTransformer(nBits=dim, **kwargs)
            else:
                raise ValueError('Invalid transformer')

            # get compound smiles from the graph
            smiles = self.get_smiles()
            # create a DataFrame from the smiles
            smiles = pd.DataFrame(smiles.items(), columns=['id', 'smiles'])
            # suppress warning from RDKit and summarize warning
            with contextlib.redirect_stderr(None):
                PandasTools.AddMoleculeColumnToFrame(smiles, smilesCol='smiles')

            # get rows id with NaN in the ROMol column
            valid_row_id = smiles.loc[~smiles['ROMol'].isna(), 'id'].tolist()
            unvalid_row_id = smiles.loc[smiles['ROMol'].isna(), 'id'].tolist()

            logging.warning(f'Successfully parse {len(smiles) - len(unvalid_row_id)} rows with valid SMILES from the MOL file!\n'
                            f'total {len(unvalid_row_id)} Invalid rows with "Unhandled" in the ROMol column:\n'
                            f' {unvalid_row_id}, removed from the final output!')

            smiles = smiles.dropna(subset=['ROMol'])
            # get the molecular vector
            mol_embeddings = transformer.transform(smiles['ROMol'])
            self.mol_embedding = dict(zip(valid_row_id, mol_embeddings))
            return self.mol_embedding


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
                clique_edges = pd.concat([clique_edges,
                                         pd.DataFrame.from_records({'entry1': member1, 'entry2': member2, 'type': 'PPrel', 'name': 'protein-complex', 'value': 'custom'}, index=[0])])

        df['entry1'] =   [group_id_dict[entry]  if entry in group_id_dict else entry for entry in df['entry1']]
        df['entry2'] =   [group_id_dict[entry]  if entry in group_id_dict else entry for entry in df['entry2']]
        # will also need add inter-interaction between the group members
        df = df.explode('entry1').explode('entry2')
        df = pd.concat([df, clique_edges])
        return df

    def process_edge_direction(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Process the direction of the edges
        '''


        # for type=PPrel, name=protein-complex, set the direction to undirected
        df.loc[(df['type'] == 'PPrel') & (df['name'] == 'protein-complex'), 'direction'] = 'undirected'
        # for type=PPrel or GErel, name="state change" or "binding/association", "dissociation", "missing interaction",
        # set the direction to undirected otherwise set the direction to directed
        undirected_PPrel = ['state change', 'binding/association', 'dissociation', 'missing interaction']
        df.loc[(df['type'].isin(['PPrel', 'GErel', 'PCrel'])) & (df['name'].isin(undirected_PPrel)), 'direction'] = 'undirected'
        # for type=PCrel and name="protein-protein expansion", set the direction to undirected
        df.loc[(df['type'] == 'PCrel') & (df['name'] == 'protein-protein expression'), 'direction'] = 'undirected'
        # for type=PCrel and name="enzyme-enzyme expansion", set the direction to ?
        df['direction'] = df['direction'].fillna('directed')


    def get_reaction_dict(self):
        '''
        Get the reaction from the graph
        '''
        reactions = [(reaction.get('id'), reaction.get('name'), reaction.get('type')) for reaction in self.root.findall('entry')]
        reactions = pd.DataFrame(reactions, columns=['id', 'name', 'type'])
        # names are separated by space, make it as a list


    def get_smiles(self) -> dict:
        if self.mol_smiles is None:
            smiles = []
            compounds = get_unique_compound_values(self.interaction)
            for compound in compounds:
                mol_file_string = fetch_mol_file_string(compound)
                smiles.append(get_smiles(mol_file_string))
                self.mol_smiles = dict(zip(compounds, smiles))
        return self.mol_smiles

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

    def get_protein_embedding(self, embedding_file):
        '''
        Get the protein embeddings
        '''
        if self.protein_embedding is None:
             self.protein_embedding = load_embedding_from_h5(embedding_file)
        return self.protein_embedding


    def get_protein_features(self):
        '''
        Get the protein features
        '''
        if self.id_converter is None or self.id_converter.target != 'uniprot':
            raise ValueError('The target must be set to uniprot')

        proteins = get_unique_proteins(self.interaction)
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



    def to_pyg_data(self):
        '''
        Convert the graph to PyG data
        '''
        pass



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