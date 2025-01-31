import os

import h5py
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import InMemoryDataset, HeteroData

class MetaboliteProteinInteraction(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ['9606.protein.links.v12.0_translated.csv',
                '9606.protein_chemical.links.v5.0.filtered.tsv',
                # 'chemical_chemical.links.detailed.v5.0.tsv.'
                'protein_embeddings.h5',
                'compound_embeddings_PCA.h5']

    @property
    def processed_file_names(self):
        return ['data.pt']

    def download(self):
        raise NotImplementedError('Data download not supported. '
                                  'Please make sure the following files are in the raw_dir: ', self.raw_file_names)

    def process(self):
        data = self.process_data(self.raw_dir, self.raw_file_names)
        data_list = [data]
        if self.pre_filter is not None:
            data_list = [data for data in data_list if self.pre_filter(data)]

        if self.pre_transform is not None:
            data_list = [self.pre_transform(data) for data in data_list]

        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

    def process_data(self, raw_dir, raw_file_names):
        # if name contain tsv, read as tsv, else read as csv
        if 'tsv' in raw_file_names[0]:
            ppi = pd.read_csv(os.path.join(raw_dir, raw_file_names[1]), sep='\t', usecols=[0, 1])
        else:
            ppi = pd.read_csv(os.path.join(raw_dir, raw_file_names[0]), usecols=[0, 1])
        ppi.columns = ['source', 'target']
        proteins = set(ppi['source']).union(ppi['target'])

        if 'tsv' in raw_file_names[1]:
            mpi = pd.read_csv(os.path.join(raw_dir, raw_file_names[1]), sep='\t', usecols=[0, 1])
        else:
            mpi = pd.read_csv(os.path.join(raw_dir, raw_file_names[1]), usecols=[0, 1])
        mpi.columns = ['source', 'target']
        compounds = set(mpi['source'])

        protein_embeddings = load_embedding_from_h5(os.path.join(raw_dir, raw_file_names[2]))
        compound_embeddings = load_embedding_from_h5(os.path.join(raw_dir, raw_file_names[3]))

        proteins = list(proteins.intersection(protein_embeddings.keys()))
        compounds = list(compounds.intersection(compound_embeddings.keys()))

        data = HeteroData()

        # Add nodes and their features
        data['gene'].x = torch.tensor([protein_embeddings[node] for node in proteins if node in protein_embeddings], dtype=torch.float)
        data['compound'].x = torch.tensor([compound_embeddings[node] for node in compounds if node in compound_embeddings], dtype=torch.float)

        # Add edges
        protein_indices = {node: i for i, node in enumerate(proteins)}
        compound_indices = {node: i for i, node in enumerate(compounds)}

        # create edge index in tensor form [2, num_edges] i.e. [[source], [target]]
        edge_index_protein = torch.tensor([[protein_indices[row['source']], protein_indices[row['target']]] for _, row in ppi.iterrows() if row['source'] in protein_indices and row['target'] in protein_indices], dtype=torch.long).t().contiguous()
        edge_index_compound = torch.tensor([[compound_indices[row['source']], protein_indices[row['target']]] for _, row in mpi.iterrows() if row['source'] in compound_indices and row['target'] in protein_indices], dtype=torch.long).t().contiguous()

        data['gene', 'to', 'gene'].edge_index = edge_index_protein
        data['compound', 'to', 'gene'].edge_index = edge_index_compound

        mapping = {'gene': protein_indices, 'compound': compound_indices}
        # save the mapping to root
        torch.save(mapping, os.path.join(self.root, 'mapping.pt'))

        return data

def load_embedding_from_h5(file_path):
    with h5py.File(file_path, 'r') as h5file:
        embeddings = {key: value[()] for key, value in h5file.items()}
    return embeddings