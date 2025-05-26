# KAMPING (KEGG Automated Metabolite Protein Interaction Network for Graph-model)

KAMPING is a powerful Python library for analyzing and modeling KEGG pathway data using graph neural networks. It provides tools for downloading, processing, and analyzing KEGG pathway data, with a focus on metabolite-protein interactions.

## Features

- Download and parse KEGG KGML files
- Create and manipulate KEGG pathway graphs
- Support for different graph types:
  - Mixed graphs (genes and metabolites)
  - Gene-only graphs
  - Metabolite-only graphs
- Graph neural network models for pathway analysis
- Protein and metabolite embedding processing
- Heterogeneous graph neural network support

## Installation

KAMPING can be installed using Poetry:

```bash
# Clone the repository
git clone https://github.com/FDUguchunhui/KAMPING.git
cd KAMPING

# Install dependencies using Poetry
poetry install

# manually insall torch
pip install torch
```

## Basic Usage

### 1. Downloading KEGG KGML Files

```python
import kamping

# Download KGML files for a specific species (e.g., human - 'hsa')
kamping.kgml('hsa', out_dir='data/kgml_hsa', verbose=True)
```

### 2. Creating a KEGG Graph

```python
# Create a mixed graph from a KGML file
mixed_graph = kamping.KeggGraph(
    'data/kgml_hsa/hsa00010.xml',
    type='mixed',
    gene_group_as_interaction=False,
    multi_substrate_as_interaction=False,
    auto_correction='fix',
    directed=True
)

# Access graph nodes and edges
print(mixed_graph.nodes)
print(mixed_graph.edges)
```

### 3. Working with Different Graph Types

```python
# Create a gene-only graph
gene_graph = kamping.KeggGraph('data/kgml_hsa/hsa00010.xml', type='gene')

# Create a metabolite-only graph
metabolite_graph = kamping.KeggGraph('data/kgml_hsa/hsa00010.xml', type='metabolite')
```

## Tutorials

The project includes several tutorials to help you get started:

1. Introduction and Basic Usage
2. Homogeneous Graph Neural Network Model
3. Protein Embedding Processing
4. Heterogeneous Graph Neural Network
5. Metabolite Embedding Processing

## Dependencies

- Python 3.10+
- NetworkX
- PyTorch
- PyTorch Geometric
- RDKit
- scikit-mol
- pandas
- numpy
- matplotlib
- and more (see pyproject.toml for complete list)

## License

This project is licensed under the terms of the included LICENSE file.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use KAMPING in your research, please cite:

[Citation information to be added]

## Contact

For questions and support, please open an issue on the GitHub repository or contact Chunhui Gu (cgu3@mdanderson.org). 