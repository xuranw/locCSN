import os
import numpy as np
import time, math
import pandas as pd

from scipy.sparse import csr_matrix, find
from scipy.stats import norm
from joblib import Parallel, delayed

def which(x):
    n = len(x)
    m = sum(x)
    y = np.zeros(m)
    j = -1
    for i in range(0, n):
        if x[i]:
            j = j+1
            y[j] = i
    y = y.astype(int)
    return y
    
def csntoflat(csn):
    if type(csn) is list:
        n2 = len(csn)
        n1 = csn[0].shape[0]
        csn_flat = np.zeros((int(n1*(n1-1)/2), n2))
        k = 0
        for i in range(0, n1-1):
            for j in range(i+1, n1):
                csn_flat[k, :] = np.asarray([item[i,j] for item in csn_loc_mat])
                k = k + 1
    else: 
        (n1, n11, n2) = csn.shape
        if n1 != n11: 
            print('dimension not match!')
            return
        csn_flat = np.zeros((int(n1*(n1-1)/2), n2))
        k = 0
        for i in range(0, n1-1):
            for j in range(i+1, n1):
                csn_flat[k, :] = csn[i, j, :]
                k = k + 1
    return csn_flat
    
def sparsetoid(sp_mtx, x = 0, y = 0):
    (id_x, id_y) = sp_mtx.nonzero()
    if len(id_x)*len(id_y) == 0:
        v = np.array([])
    else:
        v = np.asarray(sp_mtx[id_x, id_y])[0]
    d = {'x_id': id_x+x, 'y_id': id_y+y, 'value': v}
    df = pd.DataFrame(data=d)
    return df  
    
def id_concat(pair):
    frames = [pair[0], pair[1]]
    result = pd.concat(frames)
    return result
    
def idtosparse(df, G1=None, G2=None):
    if G1 is None:
        G1 = max(df.x_id) + 1
    if G2 is None:
        G2 = max(df.y_id) + 1
    sp_mtx = csr_matrix((df.value, (df.x_id, df.y_id)), shape=(G1, G2))
    return sp_mtx
  
def valuetosparse(value, I, J, G1 = None, G2 = None):
    if G1 is None:
        G1 = max(df.x_id) + 1
    if G2 is None:
        G2 = max(df.y_id) + 1
    I = I[value != 0]
    J = J[value != 0]
    value = value[value != 0]
    sp_mtx = csr_matrix((value, (I, J)), shape = (G1, G2))
    return sp_mtx
