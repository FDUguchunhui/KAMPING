import sys
from typing import List, Dict

import h5py
import numpy as np
import pandas as pd
import os.path as osp
import os

import torch
from overrides import overrides

import warnings
from torch_geometric.data import Data




class KeggGraphDataset(Data):


    def __init__(self, transform=None, pre_transform=None, pre_filter=None,):

        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = torch.load(self.processed_paths[0])
