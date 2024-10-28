import copy
import logging
import os
from typing import Union, Any, List, Dict, Literal

import h5py
import networkx as nx
import numpy as np
import pandas as pd
import requests
import scikit_mol.fingerprints
import torch
import torch_geometric as pyg
import typer
from rdkit import Chem as Chem

from kamping.parser.utils import load_embedding_from_h5, get_unique_compound_values
from kamping.data.convert import from_hetero_networkx

def load_node_h5(file_path: str):

    with h5py.File(file_path, 'r') as file:
        nodes = file.keys()
        mapping = {index: i for i, index in enumerate(nodes)}
        embedding = [torch.tensor(value) for value in file.values()]
        embedding = torch.stack(embedding)

    return embedding, mapping


def convert_to_pyg(graph, embeddings, unmatch_embeddings=None, verbose=True) -> pyg.data.Data:
    '''
    Convert the graph to PyG data
    '''
    logger = logging.getLogger()
    if verbose:
        logger.setLevel(logging.INFO)

    if not embeddings:
        raise ValueError('empty embeddings information given!')
    if graph.type == 'gene' or graph.type == 'metabolite':
        # convert the graph to networkx
        G = graph.to_networkx()
        node_without_embedding = [node for node in G.nodes if node not in embeddings.keys()]
        if node_without_embedding:
            G.remove_nodes_from(node_without_embedding)
        # add embedding to the node attributes
        nx.set_node_attributes(G, embeddings, name='embeddings')
        if graph.type == 'gene':
            node_names = dict(zip(graph.genes, graph.genes))
        else:
            node_names = dict(zip(graph.compounds, graph.compounds))
        nx.set_node_attributes(G, node_names, name='node_name')
        # convert the graph to PyG data
        if nx.is_empty(G):
            raise ValueError(f'{graph.name} is empty using giving setting!')
        data =  pyg.utils.from_networkx(G, group_node_attrs=['embeddings'])
    else:
        raise ValueError('Graph type must be "gene" or "metabolite"!')
    # logging succesful convert
    logging.info(f'{graph.name} converted to torch_geometric successfully')
    return data



def convert_to_hetero_pyg(graph, protein_embeddings, mol_embeddings):
    if graph.type != 'mixed':
        raise ValueError('graph type must be mixed!')
    G = graph.to_networkx()
    # combine protein_embeddings and mol_embeddings
    embeddings = {**protein_embeddings, **mol_embeddings}
    nx.set_node_attributes(G, embeddings, name='embeddings')
    data = from_hetero_networkx(G, node_type_attribute='node_type',
                                group_node_attrs=['embeddings'])
    return data



def get_uniprot_protein_embeddings(graphs: Union[Any, List[Any]], embedding_file) -> dict:
    '''
    Get the protein embeddings
    '''
    # if pass a single item
    if not isinstance(graphs, list):
        graphs = [graphs]

    protein_id_types = [graph.gene_id_type for graph in graphs]
    if not all([protein_id_type == 'uniprot' for protein_id_type in protein_id_types]):
        raise ValueError('Protein ID is not uniprot! Please use converter to convert protein ID into uniprot')

    embeddings = load_embedding_from_h5(embedding_file)
    # combine all proteins from the graph in the list into one
    unique_proteins = set()
    for graph in graphs:
        unique_proteins.update(graph.genes)

    # get the protein embeddings
    protein_embeddings = {protein: embeddings[protein.removeprefix('up:')] for protein in unique_proteins if protein.removeprefix('up:') in embeddings.keys()}
    return protein_embeddings


def combine_graphs(pyg_data: list):
    '''
    Given a list of pyg graph data, combine them into one large graph data
    works for both homogenous and heterogenous graph
    '''
    # use update instance function
    data = copy.copy(pyg_data[0])
    for i in range(1, len(pyg_data)):
        data.update(pyg_data[i])

    return data


def get_smiles(interaction:pd.DataFrame) -> dict:
    smiles = []
    compounds = get_unique_compound_values(interaction)
    for compound in compounds:
        # this is inefficient no need to convert back and forth
        mol_file_string = fetch_mol_file_string(compound)
        smiles.append( Chem.MolToSmiles(Chem.MolFromMolBlock(mol_file_string)))
    return dict(zip(compounds, smiles))

def get_kegg_mol(graphs: Union[Any, list[Any]]) -> pd.DataFrame:
    '''Return a dataframe with first column compound id and second column contain mol class instance'''

    # check the compound ID type for all graphs is "kegg"
    compound_id_types = [graph.compound_id_type for graph in graphs]
    if not all([compound_id_type == 'kegg' for compound_id_type in compound_id_types]):
        raise ValueError('Compound ID is not kegg! Consider recreate graphs without compound ID conversion')

    # combine all proteins from the graph in the list into one
    unique_compounds = set()
    if not isinstance(graphs, list):
        unique_compounds = graphs.compounds
    else:
        for graph in graphs:
            unique_compounds.update(graph.compounds)

    # get the molecule from the unique compounds and return a DataFrame
    mols = get_molecule(unique_compounds, mol_column='ROMol')
    return mols

def get_mol_embeddings_from_dataframe(mols, transformer, dim=1024, **kwargs) -> dict[str, np.array]:
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

    # get rows id with NaN in the ROMol column
    valid_row_id = mols.loc[~mols['ROMol'].isna(), 'id'].tolist()
    unvalid_row_id = mols.loc[mols['ROMol'].isna(), 'id'].tolist()

    logging.warning(f'''Successfully parse {len(mols) - len(unvalid_row_id)} rows with valid SMILES from the MOL file!\n'
                    total {len(unvalid_row_id)} Invalid rows with "None" in the ROMol column''')
    if not unvalid_row_id:
        logging.warning(f' {unvalid_row_id}, removed from the final output!')

    smiles = mols.dropna(subset=['ROMol'])
    # get the molecular vector
    mol_embeddings = transformer.transform(smiles['ROMol'])
    mol_embeddings = dict(zip(valid_row_id, mol_embeddings))
    return mol_embeddings


def get_mol_embeddings(graphs: Union[Any, list[Any]], transformer, dim=1024, **kwargs) -> dict[str, np.array]:
    '''
    Get the molecular embeddings
    '''

    # get dataframe with a column has mol instance
    mols = get_kegg_mol(graphs)
    return get_mol_embeddings_from_dataframe(mols, transformer, **kwargs)


def get_molecule(unique_compounds, mol_column='ROMol', verbose=True) -> pd.DataFrame:
    # if save_path is not None create the path
    mols = []
    for compound in unique_compounds:
        # this is inefficient no need to convert back and forth
        mol_file_string = fetch_mol_file_string(compound)
        if verbose:
            typer.echo(f'Fetching {compound} from KEGG')
        mols.append(Chem.MolFromMolBlock(mol_file_string))
    # create a pd.DataFrame
    data = pd.DataFrame({'id': list(unique_compounds), mol_column: mols})
    return data


def fetch_mol_file_string(compound_id):
    '''
    Fetch the MOL file string for a given compound ID.
    '''
    url = f"https://rest.kegg.jp/get/{compound_id}/mol"
    molFile = requests.get(url)
    return(str(molFile.text))
