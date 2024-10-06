import logging
import os
import pickle

import h5py
import pandas as pd
from bio_embeddings.embed import ProtTransT5XLU50Embedder

logging.basicConfig(level=logging.INFO)

def protein_embedding(model_directory: str, sequence_file: str, output_file: str,
                      half_model: bool = True) -> None:
    '''
    Embed protein sequences using ProtTrans T5-XL-U50 model
    '''
    # Load the ProtTrans T5-XL-U50

    embedder =ProtTransT5XLU50Embedder(model_directory=model_directory, half_model=half_model)
    # embeddings = list(embedder.embed_many(["MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"]))

    df = pd.read_csv('data/sequences.tsv', sep='\t', names=['id', 'sequence'], header=0)
    # remove if sequence is empty
    df = df[df['sequence'].notna()]

    embeddings = list(embedder.embed_many(df['sequence'].tolist()))

    # protein_seqs = pd.read_csv('data/sequences.tsv', sep='\t')['sequence'].tolist()
    # embedding = list(embedder.embed_many(protein_seqs))
    embeddings = [embedder.reduce_per_protein(e) for e in embeddings]

    # save protein id and embedding as a dictionary
    embeddings = dict(zip(df['id'], embeddings))

    # create the output directory if it does not exist
    os.makedirs(os.path.dirname('data/embedding/protein_embedding.h5'), exist_ok=True)

    # save the embedding to as h5 file
    with h5py.File(output_file, 'w') as h5file:
        for key, value in embeddings.items():
            h5file.create_dataset(key, data=value)
        logging.info(f"Embedding saved to {output_file}")


if __name__ == '__main__':
    protein_embedding('/Users/cgu3/embedding_model/prottrans_t5_xl_u50', 'data/sequences.tsv', 'data/embedding/protein_embedding.h5')