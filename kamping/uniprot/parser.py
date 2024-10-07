import logging

import pandas as pd
from typing_extensions import deprecated
from unipressed import UniprotkbClient

import kamping.utils

logging.basicConfig(level = logging.INFO)
import os
from typing import List

import requests

import kamping.uniprot.uniprot


def get_protein_seq(accession_number: str) -> str:
    '''
    Retrieve a protein sequence from UniProtKB using accession number

    Parameters:
    accession_number (str): UniProtKB accession number

    Returns:
    str: Protein sequence
    '''

    url = f"https://rest.uniprot.org/uniprotkb/search?query={accession_number}&fields=sequence&format=tsv"
    response = requests.get(url)

    if response.status_code == 200:
        tsv_data = response.text
        # Extract the sequence from the TSV format
        sequence = tsv_data.split('\n')[1]
        return sequence
    else:
        raise ValueError(f"Failed to retrieve data for accession number {accession_number}")

def get_protein_seqs(accession_numbers: List[str]) -> List[str]:
    '''
    Retrieve multiple protein sequences from UniProtKB using accession numbers

    Parameters:
    accession_numbers (List[str]): List of UniProtKB accession numbers

    Returns:
    List[str]: List of protein sequences
    '''
    sequences = []
    for accession_number in accession_numbers:
        logging.info(f"Retrieving sequence for {accession_number}")
        sequence = get_protein_seq(accession_number)
        sequences.append(sequence)
    return sequences

def get_protein_seqs_from_files(folder, output_path, batch:int=500) -> None:
    '''
    Retrieve protein sequences from files in a folder and save to a file

    This could be too slow
    '''
    # if folder is empty raise error
    if not os.listdir(folder):
        raise ValueError(f"Folder {folder} is empty")
    df = kamping.utils.read_all_tsv_files(folder)
    proteins = kamping.uniprot.uniprot.get_unique_proteins(df, 'up')

    # todo: a naive for batch implementation
    # check Uniprot pagination for a better implementation
    all_seqs = []
    for i in range(0, len(proteins), batch):
        # progress  bar
        logging.info(f"Fetching batch {i//batch + 1} of {len(proteins)//batch}")
        batch_proteins = proteins[i:i + batch]
        res = UniprotkbClient.fetch_many(batch_proteins)
        seqs = [entry['sequence']['value'] for entry in res]
        all_seqs.extend(seqs)


    # save as a tsv file
    df = pd.DataFrame({'entry': proteins, 'sequence': all_seqs})

    # create the output directory if it does not exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    get_protein_seqs_from_files(folder='data/converted', output_path='data/sequences.tsv')