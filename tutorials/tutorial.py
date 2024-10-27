import os
import kamping.parser.network as network
# get list of all files path in data/kgml_hsa
directory = 'data/kgml_hsa'
file_paths = [os.path.join(directory, file) for file in os.listdir(directory)]
# apply network.KeggGraph() to each file
graphs = [network.KeggGraph(file, type='metabolite') for file in file_paths]