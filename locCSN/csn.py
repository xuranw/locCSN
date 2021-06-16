import os
import numpy as np
import time, math
import pandas as pd

from scipy.sparse import csr_matrix, find
from scipy.stats import norm
from joblib import Parallel, delayed

def csn(data_full, g_mtx = None, wd_q = 0.1, dev = True, md = 1, iteration = False, fuzzy = False, ncore = 4):
    (n1, n2) = data_full.shape
    eps = np.finfo(float).eps
    if g_mtx is None:
        g_mtx = np.ones((n1, n1)) - np.tri(n1)
        zero_id = np.where(data_full.sum(axis = 1)==0)
        g_mtx[zero_id, :] = 0
        g_mtx[:, zero_id] = 0
    #csn = [[[0 for col in range(n1)]for row in range(n1)] for x in range(n2)]
    #csn = np.zeros((n1, n1, n2))
    
    (I, J, S) = find(g_mtx)
    L = len(I)
    print(*[L , 'pairs need calculation'])
    csn_mat = np.zeros((L, n2))
    
    if dev:
        if fuzzy:
            selected_gene = np.unique([I, J])
            grid_gene = np.zeros((n1, n2))
            gene_first_zero = np.zeros(n1)
            gene_mid = np.zeros((n1, 20))
            
            for s in selected_gene:
                gene_temp = data_full[s, :]
                max_gene = max(gene_temp[gene_temp!=0])
                min_gene = min(gene_temp[gene_temp!=0])
                range_temp = max_gene - min_gene
                gene_cut = np.arange(min_gene, max_gene, range_temp/20)
                gene_mid[s,:] = gene_cut[0:20] + range_temp/40
                if sum(gene_temp == 0) > 1:
                    gene_cut = np.insert(gene_cut, 0, 0)
                    gene_first_zero[s] = 1
                grid_gene[s,:] = np.digitize(gene_temp, gene_cut)
            
            def inner_fuzzy_fun(m):
                i = I[m]
                j = J[m]
                gene1 = data_full[i, :]
                gene2 = data_full[j, :]
                data = data_full[[i, j], :]
                grid1 = grid_gene[i, :]
                grid2 = grid_gene[j, :]
                grid_mix = grid1*30 + grid2
                u_grid = np.unique(grid_mix)
                n_grid = np.zeros(len(u_grid))
                for s in range(0, len(u_grid)):
                    n_grid[s] = sum(grid_mix == u_grid[s])
                
                u_grid_mix = np.vstack((np.floor((u_grid - 1)/30), (u_grid - 1)%30 + 1, u_grid, n_grid))
                if gene_first_zero[i] == 1:
                    u_grid_mix = u_grid_mix[:, u_grid_mix[0, :]!=1]
                if gene_first_zero[j] == 1:
                    u_grid_mix = u_grid_mix[:, u_grid_mix[1, :]!=1]
                cell_id = np.zeros(u_grid_mix.shape[1])
                cell_id_full = []
                for t in range(0, u_grid_mix.shape[1]):
                    cell_id_full.append(np.where(grid_mix == u_grid_mix[2, t])[0])
                    cell_id[t] = np.random.choice(cell_id_full[t], 1)
                cell_id = cell_id.astype(int)
                (upper, lower) = upperlower_dev(gene1, gene2, wd_q, md, iteration, cell_id)
                csn_temp = np.zeros(n2)
                for t in range(0, len(cell_id)):
                    k = cell_id[t]
                    B = np.zeros(2, n2)
                    for l in range(0, n2):
                        B[:, l] = (data[:,l] <= upper[:, k]) & (data[:, l] >= lower[:, k]) & (data[:, k] > 0)
                    a = B.sum(axis = 1)
                    a = np.reshape(a, (2, 1))
                    temp = (B@B.T*n2 - a@a.T)/np.sqrt((a@a.T)*((n2-a)@(n2-a).T)/(n2-1)+eps)
                    csn_temp[cell_id_full[t]] = temp[0, 1]
                
                return csn_temp
            csn_arr_temp = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fuzzy_fun)(m) for m in range(0, L)))
            csn = [valuetosparse(v, I, J, n1, n1) for v in list(csn_arr_temp.T)]
        else:
            def inner_fun(m):
                i = I[m]
                j = J[m]
                gene1 = data_full[i,:]
                gene2 = data_full[j,:]
                data = data_full[[i,j],:]
                csn_temp = np.zeros(n2)
                (upper, lower) = upperlower_dev(gene1, gene2, boxsize = wd_q, md = md, iteration = iteration)
                for k in range(0, n2):
                    if gene1[k]*gene2[k] > 0:
                        B = np.zeros((2, n2))
                        for l in range(0, n2):
                            B[:, l] = (data[:,l] <= upper[:, k]) & (data[:, l] >= lower[:, k]) & (data[:, k] > 0)
                        a = B.sum(axis = 1)
                        a = np.reshape(a, (2, 1))
                        temp = (B@B.T*n2 - a@a.T)/np.sqrt((a@a.T)*((n2-a)@(n2-a).T)/(n2-1)+eps)
                        csn_temp[k] = temp[0, 1]
                return csn_temp
            csn_arr_temp = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fun)(m) for m in range(0, L)))
            csn = [valuetosparse(v, I, J, n1, n1) for v in list(csn_arr_temp.T)]
            
    else: 
        (upper, lower) = upperlower(data_full, boxsize = wd_q)
        csn = []
        for k in range(0, n2):
            B = np.zeros((n1, n2))
            for j in range(0, n2):
                B[:, j] = (data_full[:, j] <= upper[:, k]) & (data_full[:, j] >= lower[:, k]) & (data_full[:, k] > 0)
            a = B.sum(axis = 1)
            a = np.reshape(a, (n1, 1))
            temp = (B@B.T*n2 - a@a.T)/np.sqrt((a@a.T)*((n2-a)@(n2-a).T)/(n2-1)+eps)
            np.fill_diagonal(temp, 0)
            csn.append(csr_matrix(temp))
    return csn
    
def upperlower_dev(gene1, gene2, boxsize = 0.1, md = 1, iteration = False, cell_id = None):
    if len(gene1) != len(gene2):
        return
    n1 = 2
    n2 = len(gene1)
    data = np.append([gene1], [gene2], axis= 0)
    if cell_id is None:
        cell_id = range(0, n2)
    (up_q, low_q) = upperlower(data, boxsize)
    upper = np.zeros((n1, n2))
    lower = np.zeros((n1, n2))
    
    if iteration:
        maxiter = 10000
        for k in cell_id:
            if gene1[k] * gene2[k] > 0:
                d2_0 = md * gene2[(gene1 <= up_q[0,k]) & (gene1 >= low_q[0, k])].std()
                d1_0 = md * gene1[(gene2 <= up_q[1,k]) & (gene2 >= low_q[1, k])].std()
                d2_1 = md * gene2[(gene1 <= gene1[k] + d1_0) & (gene1 >= gene1[k] - d1_0)].std()
                d1_1 = md * gene1[(gene2 <= gene2[k] + d2_0) & (gene2 >= gene2[k] - d2_0)].std()
                count = 0
                while (math.sqrt(pow(d2_0-d2_1, 2)+pow(d1_0-d1_1, 2)) < pow(10, -5)) & count < maxiter:
                    d2_0 = d2_1
                    d1_0 = d1_1
                    d2_1 = md * gene2[(gene1 <= gene1[k] + d1_0) & (gene1 >= gene1[k] - d1_0)].std()
                    d1_1 = md * gene1[(gene2 <= gene2[k] + d2_0) & (gene2 >= gene2[k] - d2_0)].std()
                    count = count + 1
                if count >= 10000:
                    print('Iteration at cell' , k, ' exceeds ', maxiter)
                    return
                upper[0, k] = gene1[k] + d1_1
                upper[1, k] = gene2[k] + d2_1
                lower[0, k] = gene1[k] - d1_1
                lower[1, k] = gene2[k] - d2_1
    else:
        for k in cell_id:
            if gene1[k] * gene2[k] > 0:
                d2 = md * gene2[(gene1 <= up_q[0,k]) & (gene1 >= low_q[0, k])].std()
                d1 = md * gene1[(gene2 <= up_q[1,k]) & (gene2 >= low_q[1, k])].std()
                upper[0, k] = gene1[k] + d1
                upper[1, k] = gene2[k] + d2
                lower[0, k] = gene1[k] - d1
                lower[1, k] = gene2[k] - d2
    return (upper, lower)
    
def upperlower(data, boxsize = 0.1):
    (n1, n2) = data.shape # n1 gene; n2 cells
    upper = np.zeros((n1, n2))
    lower = np.zeros((n1, n2))
    
    for i in range(0, n1):
        s1 = sorted(data[i,:])
        s2 = data[i,:].argsort()
        #s1.append(0)
        h = round(boxsize/2 * n2)
        k = 0
        while k < n2:
            s = 0    
            #while (k + s + 1 < n2) & (s1[k + s + 1] == s1[k]):
            while k+s+1 < n2:
                if s1[k+s+1] == s1[k]:
                    s = s + 1
                else:
                    break
            
            if s >= h:
                upper[i, s2[k:k + s + 1]] = data[i, s2[k]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[k]]
            else:
                upper[i, s2[k:k + s + 1]] = data[i, s2[min(n2 - 1, k + s + h)]]
                lower[i, s2[k:k + s + 1]] = data[i, s2[max(0, k - h)]]
    
            k = k + s + 1
    return (upper, lower)     
    
    
def upperlower_soft(data_full, soft_c, wd_q = 0.1):
    (n2, K) = soft_c.shape
    (n1, n2) = data_full.shape
    F_c = soft_c/sum(soft_c)
    upper = [np.zeros((n1, n2)), np.zeros((n1, n2))]
    lower = [np.zeros((n1, n2)), np.zeros((n1, n2))]
    for cl in range(0, K):
        fc = F_c[:, cl]
        for i in range(0, n1):
            s1 = sorted(data_full[i,:])
            s2 = data_full[i,:].argsort()
            n3 = s1.count(0)
            k = 0
            while k < n2:
                s = 0
                while k+s+1 < n2:
                    if s1[k+s+1] == s1[k]:
                        s = s + 1
                    else:
                        break
                if sum(fc[s2[k:k+s+1]]) >= wd_q/2:
                    upper[cl][i, s2[k:k+s+1]] = data_full[i, s2[k]]
                    lower[cl][i, s2[k:k+s+1]] = data_full[i, s2[k]]
                else:
                    h = 1
                    while (h+k+s < n2) & (sum(fc[s2[k:k+s+h+1]]) < wd_q/2):
                        h = h+1
                    upper[cl][i, s2[k:k+s+1]] = data_full[i, s2[min(n2-1, k+h+s)]]
                    h = 1
                    while (k-h>=0) & (sum(fc[s2[k-h:k+1]]) < wd_q/2):
                        h = h+1
                    lower[cl][i, s2[k:k+s+1]] = data_full[i, s2[max(n3*(n3>h), k-h)]]
                k = k + s + 1
            print('soft cluster', cl+1, 'gene', i , 'is done!')
    return(upper, lower)
                
def csn_soft_dev(data_full, soft_c, upper = None, lower = None, wd_q = 0.1, md = 1, iteration = False, maxiter = 10000):
    if upper is None or lower is None:
        (upper, lower) = upperlower_soft(data_full, soft_c, wd_q)
    K = soft_c.shape[1]
    (n1, n2) = data_full.shape
    csn = [np.zeros((n1, n1, n2)), np.zeros((n1, n1, n2))]
    for cl in range(0, K):
        n_cl = sum(soft_c[:, cl])
        soft_cl = soft_c[:, cl]
        for i in range(0, n1):
            for j in range(i+1, n1):
                nz_index = which(data_full[i,:]*data_full[j,:]*soft_cl > 0)
                for k in nz_index:
                    btw_i = (data_full[i,:] <= upper[cl][i, k]) & (data_full[i, :] >= lower[cl][i, k])
                    btw_j = (data_full[j,:] <= upper[cl][j, k]) & (data_full[j, :] >= lower[cl][j, k])
                    sdj_0 = md*np.sqrt(np.cov(data_full[j,btw_i], aweights=soft_cl[btw_i]))
                    sdi_0 = md*np.sqrt(np.cov(data_full[i,btw_j], aweights=soft_cl[btw_j]))
                    if iteration:
                        btw_i = (data_full[i,:] <= data_full[i, k]+sdi_0) & (data_full[i, :] >= data_full[i, k]-sdi_0)
                        btw_j = (data_full[j,:] <= data_full[j, k]+sdj_0) & (data_full[j, :] >= data_full[j, k]-sdj_0)
                        sdj_1 = md*np.sqrt(np.cov(data_full[j,btw_i], aweights=soft_cl[btw_i]))
                        sdi_1 = md*np.sqrt(np.cov(data_full[i,btw_j], aweights=soft_cl[btw_j]))
                        count = 0
                        while ((sdi_0-sdi_1)**2 + (sdj_0-sdj_1)**2 > pow(10, -12)) & (count < maxiter) & (sdi_1*sdj_1 >0):
                            sdi_0 = sdi_1
                            sdj_0 = sdj_1
                            btw_i = (data_full[i,:] <= data_full[i, k]+sdi_0) & (data_full[i, :] >= data_full[i, k]-sdi_0)
                            btw_j = (data_full[j,:] <= data_full[j, k]+sdj_0) & (data_full[j, :] >= data_full[j, k]-sdj_0)
                            sdj_1 = md*np.sqrt(np.cov(data_full[j,btw_i], aweights=soft_cl[btw_i]))
                            sdi_1 = md*np.sqrt(np.cov(data_full[i,btw_j], aweights=soft_cl[btw_j]))
                            count = count + 1
                        if count >= maxiter:
                            print('Iteration of Cluster', cl+1, 'at gene', i+1, 'and gene', j+1, 'at cell', k+1, 'have exceed ', maxiter)
                            return
                        sdj = sdj_1
                        sdi = sdi_1
                    else:
                        sdj = sdj_0
                        sdi = sdi_0
                    
                    nx = soft_cl[(data_full[i, :]<= data_full[i,k]+sdi) & (data_full[i, :] >= data_full[i,k]-sdi)].sum()
                    ny = soft_cl[(data_full[j, :]<= data_full[j,k]+sdj) & (data_full[j, :] >= data_full[j,k]-sdj)].sum()
                    nxy = soft_cl[(data_full[i, :]<= data_full[i,k]+sdi) & (data_full[i, :] >= data_full[i,k]-sdi) & (data_full[j, :]<= data_full[j,k]+sdj) & (data_full[j, :] >= data_full[j,k]-sdj)].sum()
                    rho_xy = nxy/n_cl - (nx/n_cl)*(ny/n_cl)
                    sigma_xy = nx*ny*(n_cl-nx)*(n_cl-ny)/(n_cl**4*(n_cl-1))
                    csn[cl][i, j, k] = rho_xy/np.sqrt(sigma_xy)
                    csn[cl][j, i, k] = rho_xy/np.sqrt(sigma_xy)
        print('soft cluster', cl+1)
        
    return csn
    
def csn_comb_cluster(csn, soft_c):
    (n2, K) = soft_c.shape
    scale = np.sqrt(soft_c[:, 1]**2 + soft_c[:, 0]**2)
    n1 = csn[0].shape[0]
    csn_comb = np.zeros((n1, n1, n2))
    for k in range(0, K):
        for i in range(0, n2):
            csn_comb[:, :, i] = (csn_comb[:, :, i] + csn[k][:, :, i]*soft_c[i, k])/scale[i]
    return csn_comb
    
def csn_rec(data1, data2, g_mtx = None, wd_q = 0.1, dev = True, md = 1, iteration = False, fuzzy = False, ncore = 4):
    (G1, N) = data1.shape
    G2 = data2.shape[0]
    eps = np.finfo(float).eps
    data = [data1, data2]
    if g_mtx is None:
        g_mtx = np.ones((G1, G2))
        zero_id1 = np.where(data1.sum(axis = 1)==0)
        g_mtx[zero_id1, :] = 0
        zero_id2 = np.where(data2.sum(axis = 1)==0)
        g_mtx[:, zero_id2] = 0
    #csn = np.zeros((G1, G2, N))
    (I, J, S) = find(g_mtx)
    L = len(I)
    print(*[L , 'pairs need calculation'])
    csn_mat = np.zeros((L, N))
    if dev:
        if fuzzy:
            selected_gene = [np.unique(I), np.unique(J)]
            grid_gene = [np.zeros((G1, N)), np.zeros((G2, N))]
            gene_first_zero = [np.zeros(G1), np.zeros(G2)]
            gene_mid = [np.zeros((G1, 20)), np.zeros((G2, 20))]
            for k in range(0, 1):
                for s in selected_gene[k]:
                    gene_temp = data[k][s, :]
                    max_gene = max(gene_temp[gene_temp!=0])
                    min_gene = min(gene_temp[gene_temp!=0])
                    range_temp = max_gene - min_gene
                    gene_cut = np.arange(min_gene, max_gene, range_temp/20)
                    gene_mid[k][s,:] = gene_cut[0:20] + range_temp/40
                    if sum(gene_temp == 0) > 1:
                        gene_cut = np.insert(gene_cut, 0, 0)
                        gene_first_zero[k][s] = 1
                    grid_gene[k][s,:] = np.digitize(gene_temp, gene_cut)
            
            def inner_fuzzy_fun(m):
                i = I[m]
                j = J[m]
                gene1 = data1[i, :]
                gene2 = data2[j, :]
                data = np.append([gene1], [gene2], axis= 0)
                grid1 = grid_gene[0][i, :]
                grid2 = grid_gene[1][j, :]
                grid_mix = grid1*30 + grid2
                u_grid = np.unique(grid_mix)
                n_grid = np.zeros(len(u_grid))
                for s in range(0, len(u_grid)):
                    n_grid[s] = sum(grid_mix == u_grid[s])
                
                u_grid_mix = np.vstack((np.floor((u_grid - 1)/30), (u_grid - 1)%30 + 1, u_grid, n_grid))
                if gene_first_zero[0][i] == 1:
                    u_grid_mix = u_grid_mix[:, u_grid_mix[0, :]!=1]
                if gene_first_zero[1][j] == 1:
                    u_grid_mix = u_grid_mix[:, u_grid_mix[1, :]!=1]
                cell_id = np.zeros(u_grid_mix.shape[1])
                cell_id_full = []
                for t in range(0, u_grid_mix.shape[1]):
                    cell_id_full.append(np.where(grid_mix == u_grid_mix[2, t])[0])
                    cell_id[t] = np.random.choice(cell_id_full[t], 1)
                cell_id = cell_id.astype(int)
                (upper, lower) = upperlower_dev(gene1, gene2, wd_q, md, iteration, cell_id)
                csn_temp = np.zeros(N)
                for t in range(0, len(cell_id)):
                    k = cell_id[t]
                    B = np.zeros(2, N)
                    for l in range(0, N):
                        B[:, l] = (data[:,l] <= upper[:, k]) & (data[:, l] >= lower[:, k]) & (data[:, k] > 0)
                    a = B.sum(axis = 1)
                    a = np.reshape(a, (2, 1))
                    temp = (B@B.T*N - a@a.T)/np.sqrt((a@a.T)*((N-a)@(N-a).T)/(N-1)+eps)
                    csn_temp[cell_id_full[t]] = temp[0, 1]
                
                return csn_temp
            csn_arr_temp = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fuzzy_fun)(m) for m in range(0, L)))
            
            #for k in range(0, n2):
            #    csn[I, J, k] = csn_mat[:, k]
            #    csn[J, I, k] = csn_mat[:, k]
        else:
            def inner_fun(m):
                i = I[m]
                j = J[m]
                gene1 = data1[i, :]
                gene2 = data2[j, :]
                data = np.append([gene1], [gene2], axis= 0)
                csn_temp = np.zeros(N)
                (upper, lower) = upperlower_dev(gene1, gene2, boxsize = wd_q, md = md, iteration = iteration)
                for k in range(0, N):
                    if gene1[k]*gene2[k] > 0:
                        B = np.zeros((2, N))
                        for l in range(0, N):
                            B[:, l] = (data[:,l] <= upper[:, k]) & (data[:, l] >= lower[:, k]) & (data[:, k] > 0)
                        a = B.sum(axis = 1)
                        a = np.reshape(a, (2, 1))
                        temp = (B@B.T*N - a@a.T)/np.sqrt((a@a.T)*((N-a)@(N-a).T)/(N-1)+eps)
                        csn_temp[k] = temp[0, 1]
                return csn_temp
            csn_arr_temp = np.asarray(Parallel(n_jobs = ncore)(delayed(inner_fun)(m) for m in range(0, L)))

        csn = [valuetosparse(v, I, J, G1, G2) for v in list(csn_arr_temp.T)]
    else: 
        print('Please use csn directly')
        return
    return csn
    
def csn_block(data, M = 100, g_mtx = None, wd_q = 0.1, dev = True, md = 1, iteration = False, fuzzy = False, ncore = 4):
    (G, K) = data.shape
    n = math.ceil(G/M)
    group_n = np.zeros(G)
    for i in range(0, n):
        group_n[i*M:min((i+1)*M, G)] = i
    if g_mtx is None:
        g_mtx = np.ones((G, G)) - np.tri(G)
        zero_id = np.where(data.sum(axis = 1)==0)
        g_mtx[zero_id, :] = 0
        g_mtx[:, zero_id] = 0
    
    csn_mtx_id = []
    for k in range(0, K):
        csn_mtx_id.append(pd.DataFrame(data={'x_id': np.array([], dtype = int), 'y_id': np.array([], dtype = int), 'value': np.array([])}))

    for i in range(0, n):
        data_i = data[group_n == i, :]
        for j in range(i, n):
            if i == j:
                g_mtx_temp = g_mtx[group_n == i, group_n == i]
                csn_temp = csn(data_i, g_mtx_temp, wd_q, dev, md, iteration, fuzzy, ncore)
                csn_id_temp = [sparsetoid(item, i*M, j*M) for item in csn_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
            else:
                data_j = data[group_n == j,:]
                g_mtx_temp = g_mtx[group_n == i, group_n == j]
                csn_temp = csn_rec(data_i, data_j, g_mtx_temp, wd_q, dev, md, iteration, fuzzy, ncore)
                csn_id_temp = [sparsetoid(item, i*M, j*M) for item in csn_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
                csn_id_temp = [item[item.columns[[1, 0, 2]]].rename(columns={'y_id': 'x_id', 'x_id': 'y_id'}) for item in csn_id_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
            print('block [', i , ',', j, '] finished!')
    
    csn_mtx = [idtosparse(item, G, G) for item in csn_mtx_id]
    return csn_mtx
    
def csn_loc(data_full, knn_index, wd_q = 0.1, dev = True, md = 1, iteration = False, ncore = 4):
    (n1, n2) = data_full.shape
    (nk, nc) = knn_index.shape
    eps = np.finfo(float).eps
    if nc != n2:
        print('dimension of data and knn do not match!')
        return
    #csn_mat = np.zeros((n1, n1, n2))
    def inner_fun(k):
        #csn_temp = np.zeros((n1, n1))
        index_temp = knn_index[:, k]
        data_sub = data_full[:, index_temp-1]
        g_index = which(data_sub[:, 0] > 0).astype(int)
        L_temp = len(g_index)
        I = np.zeros(L_temp*(L_temp-1))
        J = np.zeros(L_temp*(L_temp-1))
        S = np.zeros(L_temp*(L_temp-1))
        r = 0
        for i in range(0, L_temp-1):
            for j in range(i+1, L_temp):
                gi = g_index[i]
                gj = g_index[j]
                gene1 = data_sub[gi, :]
                gene2 = data_sub[gj, :]
                data = np.append([gene1], [gene2], axis= 0)
                if dev:
                    (upper, lower) = upperlower_dev(gene1, gene2, wd_q, md, iteration, np.array([1]))
                else:
                    (upper, lower) = upperlower(data, wd_q)
                B = np.zeros((2, nk))
                for m in range(0, nk):
                    B[:, m] = (data[:,m] <= upper[:, 1])&(data[:,m] >= lower[:, 1])
                a = B.sum(axis = 1)
                a = np.reshape(a, (2, 1))
                temp = (B@B.T*nk - a@a.T)/np.sqrt((a@a.T)*((nk-a)@(nk-a).T)/(nk-1)+eps)
                I[r] = gi
                J[r] = gj
                S[r] = temp[0, 1]
                r = r+1
                I[r] = gj
                J[r] = gi
                S[r] =temp[0, 1]
                r = r + 1
        I = I.astype(int)
        J = J.astype(int)
        I = I[S != 0]
        J = J[S != 0]
        S = S[S != 0]
        csn_temp = csr_matrix((S, (I, J)), shape=(n1, n1))
        return csn_temp
    csn_mat_list = Parallel(n_jobs = ncore)(delayed(inner_fun)(k) for k in range(0, n2))
    #for k in range(0, n2):
    #    csn_mat[:, :, k] = csn_mat_list[k]
    return csn_mat_list # csn_mat
    
def csn_rec_loc(data1, data2, knn_index, wd_q = 0.1, dev = True, md = 1, iteration = False, ncore = 4):
    (G1, N) = data1.shape
    G2 = data2.shape[0]
    (nk, nc) = knn_index.shape
    if nc != N:
        print('dimension of data and knn do not match!')
        return
    #csn_mat = np.zeros((G1, G2, N))
    def inner_fun(k):
        #csn_temp = np.zeros((G1, G2))
        index_temp = knn_index[:, k]
        data1_sub = data1[:, index_temp-1]
        data2_sub = data2[:, index_temp-1]
        g1_index = which(data1_sub[:, 0] > 0)
        L1_temp = len(g1_index)
        g2_index = which(data2_sub[:, 0] > 0)
        L2_temp = len(g2_index)
        I = np.zeros(L1_temp * L2_temp)
        J = np.zeros(L1_temp * L2_temp)
        S = np.zeros(L1_temp * L2_temp)
        r = 0
        for i in range(0, L1_temp):
            for j in range(0, L2_temp):
                gi = g1_index[i]
                gj = g2_index[j]
                gene1 = data1_sub[gi, :]
                gene2 = data2_sub[gj, :]
                data = np.append([gene1], [gene2], axis= 0)
                if dev:
                    (upper, lower) = upperlower_dev(gene1, gene2, wd_q, md, iteration, np.array([1]))
                else:
                    (upper, lower) = upperlower(data, wd_q)
                B = np.zeros((2, nk))
                for m in range(0, nk):
                    B[:, m] = (data[:,m] <= upper[:, 1])&(data[:,m] >= lower[:, 1])
                a = B.sum(axis = 1)
                a = np.reshape(a, (2, 1))
                temp = (B@B.T*nk - a@a.T)/np.sqrt((a@a.T)*((nk-a)@(nk-a).T)/(nk-1)+eps)
                I[r] = gi
                J[r] = gj
                S[r] = temp[0, 1]
                r = r+1
        #csn_temp = csr_matrix((S, (I, J)), shape=(G1, G2))
        I = I.astype(int)
        J = J.astype(int)
        I = I[S != 0]
        J = J[S != 0]
        S = S[S != 0]
        #csn_temp[I, J] = S
        csn_temp = csr_matrix((S, (I, J)), shape=(G1, G2))
        return csn_temp
    csn_mat_list = Parallel(n_jobs = ncore)(delayed(inner_fun)(k) for k in range(0, N))
    #for k in range(0, N):
    #    csn_mat[:, :, k] = csn_mat_list[k]
    return csn_mat_list # csn_mat
    
def csn_block_loc(data, knn_index, M = 100, wd_q = 0.1, dev = True, md = 1, iteration = False, ncore = 4):
    (G, K) = data.shape
    n = math.ceil(G/M)
    group_n = np.zeros(G)
    for i in range(0, n):
        group_n[i*M:min((i+1)*M, G)] = i
    
    csn_mtx_id = []
    for k in range(0, K):
        csn_mtx_id.append(pd.DataFrame(data={'x_id': np.array([], dtype = int), 'y_id': np.array([], dtype = int), 'value': np.array([])}))
    for i in range(0, n):
        data_i = data[group_n == i, :]
        for j in range(i, n):
            if i == j:
                csn_temp = csn_loc(data_i, knn_index, wd_q, dev, md, iteration, ncore)
                csn_id_temp = [sparsetoid(item, i*M, j*M) for item in csn_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
            else:
                data_j = data[group_n == j,:]
                csn_temp = csn_rec_loc(data_i, data_j, knn_index, wd_q, dev, md, iteration, ncore)
                csn_id_temp = [sparsetoid(item, i*M, j*M) for item in csn_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
                csn_id_temp = [item[item.columns[[1, 0, 2]]].rename(columns={'y_id': 'x_id', 'x_id': 'y_id'}) for item in csn_id_temp]
                csn_mtx_id = [id_concat(pair) for pair in zip(csn_mtx_id, csn_id_temp)]
            print('block [', i , ',', j, '] finished!')
    csn_mtx = [idtosparse(item, G, G) for item in csn_mtx_id]
    return csn_mtx
   
