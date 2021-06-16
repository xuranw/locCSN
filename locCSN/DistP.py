import os
import numpy as np
import time, math
import pandas as pd

from scipy.sparse import csr_matrix, find
from scipy.stats import norm
from joblib import Parallel, delayed

def create_D_mat(csn_mat, ncore = 4):
    if type(csn_mat) is list:
        T = len(csn_mat)
    else:
        T = csn_mat.shape[2]
    D_mat = np.zeros((T, T))
    s_mtx = np.ones((T, T)) - np.tri(T)
    (I, J, S) = find(s_mtx)
    L = len(I)
    if type(csn_mat) is list:
        def inner_fun(m):
            i = I[m]
            j = J[m]
            a = sum(abs(csn_mat[i]-csn_mat[j])).max()
            return a
        A = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fun)(m) for m in range(0, L)))
        D_mat[I, J] = A
        D_mat[J, I] = A
    else:
        def inner_fun(m):
            i = I[m]
            j = J[m]
            a = max(sum(abs(csn_mat[:, :, i] - csn_mat[:, :, j])))
            return a
        A = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fun)(m) for m in range(0, L)))
        D_mat[I, J] = A
        D_mat[J, I] = A
    return D_mat
    
def distance_test(D_mat, change, num_perm = 1000):
    stat_rec = np.zeros(num_perm)
    n = D_mat.shape[0]
    Z = np.zeros((n, 2))
    Z[0:change, 0] = 1
    Z[change:, 1] = -1
    
    counts = sum(abs(Z))
    counts = np.reshape(counts, (1, 2))
    count_mat = counts.T@counts
    for i in range(0, num_perm):
        if i == 0:
            ind = np.asarray(list(range(0, n)))
        else:
            ind = np.random.permutation(n)
        Dp = D_mat[ind, :][:, ind]
        test_stat = Z.T@Dp@Z / (count_mat + pow(10, -6))
        stat_rec[i] = sum(sum(test_stat))
    
    pval = 1 - (stat_rec > stat_rec[0] + pow(10, -6)).mean()
    return pval
