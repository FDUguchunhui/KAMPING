import logging
import os
from typing import Union

import networkx as nx
import pandas as pd
import torch_geometric as pyg
import kamping
from kamping import from_hetero_networkx, from_networkx


def read_all_tsv_files(directory):
    '''
    Read all files in the directory as one dataframe.
    '''
    # list a files except those start with "."

    all_files = os.listdir(directory)
    df = pd.DataFrame()
    for file in all_files:
        if file.startswith('.'):
            continue
        temp_df = pd.read_csv(os.path.join(directory, file), sep='\t')
        # add a new column "source" to the DataFrame and set it to the file name
        temp_df['source'] = file.removesuffix('.tsv')
        df = pd.concat([df, temp_df])
    # reset the index of the DataFrame
    df.reset_index(drop=True, inplace=True)
    return df


def create_graphs(directory, type, ignore_file:Union[list, None]=None, verbose=True, **kwargs):
    files = os.listdir(directory)
    if ignore_file is not None:
        for file in files:
            if file in ignore_file:
                files.remove(file)
    file_paths = [os.path.join(directory, file) for file in files]
    file_paths.sort() # in-place sort
    graphs = []
    for file in file_paths:
        try:
            graph = kamping.KeggGraph(file, type=type, verbose=verbose, **kwargs)
            graphs.append(graph)
        except Exception as e:
            logging.warning(e)
    return graphs


def convert_to_pyg(graphs: list, embeddings):
    pyg_graphs = []
    for graph in graphs:
        try:
            pyg_graphs.append(kamping.convert_to_pyg(graph, embeddings=embeddings, verbose=False))
        except Exception as e:
            logging.warning(e)
    return pyg_graphs

def convert_to_single_pyg(graphs, embeddings):
    # convert all graphs into networkx graphs
    nx_graphs = [graph.to_networkx() for graph in graphs]
    # compose graphs into one graph
    G = nx.compose_all(nx_graphs)
    G.graph['name'] = 'combined'
    graph_types = list(set([graph.type for graph in graphs]))
    all_graph_types_equal_single = all(elements == graph_types[0] for elements in graph_types)
    G.graph['type'] = graph_types[0]
    all_proteins = list(set([protein for graph in graphs for protein in graph.genes]))
    all_compounds = list(set([compound for graph in graphs for compound in graph.compounds]))

    if not embeddings:
        raise ValueError('empty embeddings information given!')
    if all_graph_types_equal_single:
        # convert the graph to networkx
        node_without_embedding = [node for node in G.nodes if node not in embeddings.keys()]
        if node_without_embedding:
            G.remove_nodes_from(node_without_embedding)
        # add embedding to the node attributes
        nx.set_node_attributes(G, embeddings, name='embeddings')
        if graph_types[0] == 'gene':
            node_names = dict(zip(all_proteins, all_proteins))
        if graph_types[0] == 'metabolite':
            node_names = dict(zip(all_compounds, all_compounds))
        else:
            node_names = dict(zip(all_proteins + all_compounds, all_proteins + all_compounds))
        nx.set_node_attributes(G, node_names, name='node_name')
        # convert the graph to PyG data
        if nx.is_empty(G):
            raise ValueError(f'Combined graph empty using giving setting!')
        if G.graph['type'] == 'gene' or G.graph['type'] == 'metabolite':
            data = from_networkx(G, group_node_attrs=['embeddings'])
        else:
            data = from_hetero_networkx(G, node_type_attribute='node_type',
                                    group_node_attrs=['embeddings'])
    else:
        raise ValueError('Graph type must be "gene", "metabolite", or "mixed"!')
    # logging succesful convert
    logging.info(f'Combined graph converted to torch_geometric successfully')
    return data