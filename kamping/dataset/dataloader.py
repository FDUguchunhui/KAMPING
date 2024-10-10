import os
import os.path as osp

import h5py
import pandas as pd
import torch

from torch_geometric.data import HeteroData, download_url, extract_zip
from torch_geometric.transforms import RandomLinkSplit, ToUndirected

import kamping
import kamping.utils


def load_node_h5(file_path: str):

    with h5py.File(file_path, 'r') as file:
        nodes = file.keys()
        mapping = {index: i for i, index in enumerate(nodes)}
        embedding = [torch.tensor(value) for value in file.values()]
        embedding = torch.stack(embedding)


    return embedding, mapping

def combined_tsv(input_dir: str):
    '''
    This function combines all the CSV files in the input directory into a single DataFrame
    '''
    # if folder is empty raise error
    if not os.listdir(input_dir):
        raise ValueError(f"Folder {input_dir} is empty")
    df = kamping.utils.read_all_tsv_files(input_dir)


def load_edge_csv(file_path: str,
                  protein_mapping:dict, mol_mapping:str):

    # if file_path is a directory, combine all the files in the directory into a single dataframe
    if os.path.isdir(file_path):
        df = kamping.utils.read_all_tsv_files(file_path)
    else:
        df = pd.read_csv(file_path, sep='\t')

    # remove rows with entry1 and entry2 both start with 'cpd:' and type is PCrel
    # remove rows with entry1 and entry2 both start with 'up:' and type is PCrel
    # additional data cleaning
    df = df[~((df['entry1'].str.startswith('cpd:')) & (df['entry2'].str.startswith('cpd:')) & (df['type'] == 'PCrel'))]
    df = df[~((df['entry1'].str.startswith('up:')) & (df['entry2'].str.startswith('up:')) & (df['type'] == 'PCrel'))]
    df = df[~((df['type'] == 'PPrel') & ((df['entry2'].str.startswith('cpd')) | (df['entry1'].str.startswith('cpd'))))]
# for df_pc if entry2 is a protein, swap entry1 and entry2
    # todo: find a more elegant way to do this
    # add it to parser function
    df.loc[(df['type'] == 'PCrel') & (df['entry2'].str.startswith('up')), ['entry1', 'entry2']] = df.loc[(df['type'] == 'PCrel') & (df['entry2'].str.startswith('up')), ['entry2', 'entry1']].values


# Remove rows with entry1 or entry2 not in the mapping
    mapping_keys = list(protein_mapping.keys()) + list(mol_mapping.keys())

    # remove prefix from entry1 and entry2
    df['entry1'] = df['entry1'].str.replace(r'^[^:]+:', '', regex=True)
    df['entry2'] = df['entry2'].str.replace(r'^[^:]+:', '', regex=True)

    # remove rows with entry1 or entry2 not in the mapping
    df = df[df['entry1'].isin(mapping_keys) & df['entry2'].isin(mapping_keys)]


    # split df into two parts: one for protein-protein edges and the other for protein-compound edges
    # the PP edges will have type "PPrel"
    # the PC edges will have type "PCrel"
    df_pp = df[df['type'] == 'PPrel']
    df_pc = df[df['type'] == 'PCrel']


    # df_pc edges
    pc_src = [protein_mapping[entry] for entry in df_pc['entry1']]
    pc_dst = [mol_mapping[entry] for entry in df_pc['entry2']]
    pc_index = torch.tensor([pc_src, pc_dst])

    # df_pp edges
    pp_src = [protein_mapping[entry] for entry in df_pp['entry1']]
    pp_dst = [protein_mapping[entry] for entry in df_pp['entry2']]
    pp_index = torch.tensor([pp_src, pp_dst])

    pc_edge_attr = None
    pp_edge_attr = None

    return pc_index, pc_edge_attr, pp_index, pp_edge_attr


protein_x, protein_mapping = load_node_h5('data/embedding/protein_embedding.h5')
mol_x, mol_mapping = load_node_h5('data/embedding/mol_embedding.h5')
pc_edge_index, _, pp_edge_index, _ = load_edge_csv('data/converted', protein_mapping, mol_mapping)

data = HeteroData()
data['protein'].num_nodes = len(protein_mapping)  # Users do not have any features.
data['protein'].node_type = 'protein'
data['protein'].x = protein_x
data['metabolite'].num_nodes = len(mol_mapping)  # Movies have features.
data['metabolite'].node_type = 'metabolite'

data['metabolite'].x = mol_x
data['protein', 'interact', 'metabolite'].edge_index = pc_edge_index
data['protein', 'interact', 'metabolite'].edge_attr = None
data['protein', 'interact', 'protein'].edge_index = pp_edge_index
data['protein', 'interact', 'protein'].edge_attr = None
print(data)

data = ToUndirected()(data)
# del data['movie', 'rev_rates', 'user'].edge_label

transform = RandomLinkSplit(
    num_val=0.05,
    num_test=0.1,
    neg_sampling_ratio=0.0,
    edge_types=[('protein', 'interact', 'metabolite')],
    rev_edge_types=[('protein', 'rev_interact', 'metabolite')],
)
train_data, val_data, test_data = transform(data)
print(train_data)
print(val_data)
print(test_data)