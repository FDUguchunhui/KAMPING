import h5py
import torch


def load_node_h5(file_path: str):

    with h5py.File(file_path, 'r') as file:
        nodes = file.keys()
        mapping = {index: i for i, index in enumerate(nodes)}
        embedding = [torch.tensor(value) for value in file.values()]
        embedding = torch.stack(embedding)


    return embedding, mapping
