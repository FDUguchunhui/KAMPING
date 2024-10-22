import copy
from typing import Union, Any, List

import h5py
import networkx as nx
import torch
import torch_geometric as pyg

from kamping.parser.utils import load_embedding_from_h5
from kamping.data.convert import from_hetero_networkx

def load_node_h5(file_path: str):

    with h5py.File(file_path, 'r') as file:
        nodes = file.keys()
        mapping = {index: i for i, index in enumerate(nodes)}
        embedding = [torch.tensor(value) for value in file.values()]
        embedding = torch.stack(embedding)

    return embedding, mapping


def convert_to_pyg(graph, embeddings) -> pyg.data.Data:
    '''
    Convert the graph to PyG data
    '''
    if not embeddings:
        raise ValueError('empty embeddings information given!')
    if graph.type == 'gene' or graph.type == 'metabolite':
        # convert the graph to networkx
        G = graph.to_networkx()
        # add embedding to the node attributes
        nx.set_node_attributes(G, embeddings, name='embeddings')
        if graph.type == 'gene':
            node_names = dict(zip(graph.proteins, graph.proteins))
        else:
            node_names = dict(zip(graph.compounds, graph.compounds))
        nx.set_node_attributes(G, node_names, name='node_name')
        # convert the graph to PyG data
        data =  pyg.utils.from_networkx(G, group_node_attrs=['embeddings'])
    else:
        raise ValueError('Graph type must be "gene" or "metabolite"!')
    return data


def convert_to_hetero_pyg(graph, protein_embeddings, mol_embeddings):
    if graph.type != 'mixed':
        raise ValueError('graph type must be mixed!')
    G = graph.to_networkx()
    # combine protein_embeddings and mol_embeddings
    embeddings = {**protein_embeddings, **mol_embeddings}
    nx.set_node_attributes(G, embeddings, name='embeddings')
    data = from_hetero_networkx(G, node_type_attribute='type',
                                group_node_attrs=['embeddings'])
    return data



def get_uniprot_protein_embeddings(graphs: Union[Any, List[Any]], embedding_file) -> dict:
    '''
    Get the protein embeddings
    '''
    # if pass a single item
    if not isinstance(graphs, list):
        graphs = [graphs]

    protein_id_types = [graph.protein_id_type for graph in graphs]
    if not all([protein_id_type == 'uniprot' for protein_id_type in protein_id_types]):
        raise ValueError('Protein ID is not uniprot! Please use converter to convert protein ID into uniprot')

    embeddings = load_embedding_from_h5(embedding_file)
    # combine all proteins from the graph in the list into one
    unique_proteins = set()
    for graph in graphs:
        unique_proteins.update(graph.proteins)

    # get the protein embeddings
    protein_embeddings = {protein: embeddings[protein] for protein in unique_proteins if protein in embeddings.keys()}
    return protein_embeddings


    def return_one_graph(pyg_data: list):
        '''
        Given a list of pyg graph data, combine them into one large graph data
        works for both homogenous and heterogenous graph
        '''
        # use update instance function
        data = copy.copy(pyg_data[0])
        for i in range(1, len(pyg_data)):
            data.update(pyg_data[i])

        return data