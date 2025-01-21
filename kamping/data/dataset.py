import os

import h5py
import networkx as nx
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import InMemoryDataset

from kamping import from_hetero_networkx


class MetaboliteProteinInteraction(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ['9606.protein_chemical_translated.links.v5.0.tsv',
                '9606.protein.links.v12.0_translated.csv',
                'protein_embeddings.h5',
                'compound_embeddings.h5']

    @property
    def processed_file_names(self):
        return ['data.pt']

    def download(self):
        raise NotImplementedError('Data download not supported.'
                                  'please make sure following files are in the raw_dir: ', self.raw_file_names)

    def process(self):
        data = self.process_data(self.raw_dir, self.raw_file_names)
        data_list = [data]
        if self.pre_filter is not None:
            data_list = [data for data in data_list if self.pre_filter(data)]

        if self.pre_transform is not None:
            data_list = [self.pre_transform(data) for data in data_list]


        # combine list of data into a big data object
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

    def process_data(self, raw_dir, raw_file_names):
        ppi = pd.read_csv(os.path.join(raw_dir, raw_file_names[0]), sep='\t', usecols=[0, 1])
        ppi.columns = ['source', 'target']
        proteins = set(ppi.iloc[:, 0].tolist() + ppi.iloc[:, 1].tolist())

        mpi = pd.read_csv(os.path.join(raw_dir, raw_file_names[1]), usecols=[0, 1])
        mpi.columns = ['source', 'target']
        compounds = set(mpi.iloc[:, 0].tolist())
        # combined two dataframes
        edges = pd.concat([ppi, mpi], ignore_index=True)

        G = nx.from_pandas_edgelist(edges,
                                    source='source', target='target', create_using=nx.DiGraph())

        node_attributes = {**{protein: 'gene' for protein in proteins},
                           **{compound: 'compound' for compound in compounds}}


        # load embeddings
        protein_embeddings = load_embedding_from_h5(os.path.join(raw_dir, raw_file_names[2]))
        compound_embeddings = load_embedding_from_h5(os.path.join(raw_dir, raw_file_names[3]))

        embeddings = {**protein_embeddings, **compound_embeddings}
        node_without_embedding = [node for node in G.nodes if node not in embeddings.keys()]
        if node_without_embedding:
            G.remove_nodes_from(node_without_embedding)

        nx.set_node_attributes(G, node_attributes, name='node_type')
        nx.set_node_attributes(G, embeddings, name='embeddings')
        data, mapping = from_hetero_networkx(G, node_type_attribute='node_type',
                                             group_node_attrs=['embeddings'])
        return data


def load_embedding_from_h5(file_path):
    '''
    Load the embedding from a h5 file
        '''
    with h5py.File(file_path, 'r') as h5file:
        embeddings = {key: value[()] for key, value in h5file.items()}
    return embeddings