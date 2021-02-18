import os
import pandas as pd
import numpy as np
import math

from scipy.sparse import csr_matrix, find
import matplotlib.pyplot as plt
import time
from scipy.stats import norm

from joblib import Parallel, delayed
import time, math

# simple example of using one cell type from Chu et al. dataset
data = pd.read_csv('logChumarker.txt', sep = ' ')
data.shape

data_temp = data.values
cell_type = pd.read_csv('chutypect.txt', sep = ' ')
celltype = cell_type.x

data_full = data_temp[:, celltype == 1];
csn_stat = csn(data_full, dev = True)  # dev = True locCSN, dev = False oCSN
# csn_stat stores test statistics

# Cutoff at norm(0.99)
csn_mat = [(item > norm.ppf(0.99)).astype(int) for item in csn_stat]
plt.imshow(sum(csn_mat).toarray()/len(csn_mat))  # plot of averaged CSN in one cell type

# Dataset from Velmeshev et al. We select 100 nearest metacells for analysis
mc_temp = pd.read_csv('mc_cpm_L.txt', sep = ' ')  # Expression
gene_name = mc_temp.index.values
mc_temp = mc_temp.values
log_mc = np.log(mc_temp+1)

meta_mc = pd.read_csv('meta_mc_L.txt', sep = ' ')  # metadata

log_temp = log_mc[:, meta_mc.cluster == 'L4']              # We focus on L4 cell group
diag_temp = meta_mc.diagnosis[meta_mc.cluster == 'L4']

knn_index = pd.read_csv('mcknn100_L4.txt', sep = ' ')     # Index of metacells that are included for locCSN construction
knn_index = knn_index.values

csn_loc_mat = csn_loc(log_temp[0:10, :], knn_index)   # Demo for small subset of genes 
csn_loc_mat = csn_block_loc(log_temp, knn_index, M = 50)     # Use parallel computation for more genes 


## Note all csn_block* function are designed for parallel computation of CSNs with large number of genes.
