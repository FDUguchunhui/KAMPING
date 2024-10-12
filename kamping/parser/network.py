#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Everest Uriel Castaneda
@desc: File for obtaining gene-only pathways
"""
import contextlib
import logging
import os
from typing import Union, Any
import h5py
import numpy as np
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
        self.root = ET.parse(input_data).getroot()
        self.graph_title = self.root.get('title')

        self.conversion_dictionary = utils.entry_id_conv_dict(self.root, unique=unique)
        self.interaction = self.get_edges()
        self.mol_smiles = None
        self.mol_embedding = None
        self.protein_embedding = None

        if self.verbose:
            typer.echo(typer.style(f"Now parsing: {self.root.get('title')}...", fg=typer.colors.GREEN, bold=False))
        # Check for compounds or undefined nodes
        has_compounds_or_undefined = not self.interaction[(self.interaction['entry1'].str.startswith('cpd:')) | (self.interaction['entry2'].str.startswith('cpd:')) | (self.interaction['entry1'].str.startswith('undefined')) | (self.interaction['entry2'].str.startswith('undefined'))].empty

        # if not mixed, remove "path" entries and propagate compounds
        if self.type == 'gene-only' :
            # Remove edges with "path" entries
            self.interaction = self.interaction[(~self.interaction['entry1'].str.startswith('path')) & (~self.interaction['entry2'].str.startswith('path'))]
            if has_compounds_or_undefined:
                self.interaction = self._propagate_compounds(self.interaction)
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

        # remove row with undefined entries
        self.interaction = self.interaction[~self.interaction['entry1'].str.startswith('undefined')]
        self.interaction = self.interaction[~self.interaction['entry2'].str.startswith('undefined')]

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

        # convert compound value to kegg id if only relation.type is "compound"
        # compound is a list with one element mapped from dict
        df['value'] =df.apply(lambda row: self.conversion_dictionary.at[row['value'], 'name'] if row['name'] == 'compound' else row['value'], axis=1)

        # add another two columns to indicate type of entry1 and entry2
        df['entry1_type'] = df['entry1'].map(lambda x: self.conversion_dictionary.at[x, 'type'])
        df['entry2_type'] = df['entry2'].map(lambda x: self.conversion_dictionary.at[x, 'type'])
        df['entry1'] = df['entry1'].map(lambda x: self.conversion_dictionary.at[x, 'name'])
        df['entry2'] = df['entry2'].map(lambda x: self.conversion_dictionary.at[x, 'name'])



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
        if self.id_converter.target != 'uniprot':
            raise ValueError('The target must be set to uniprot')

        proteins = get_unique_proteins(self.interaction)
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=proteins
        )

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