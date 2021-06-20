Local Cell-Specific Network (locCSN)
=============================================

`locCSN` is a deconvolution method that utilizes cross-subject scRNA-seq to estimate cell type proportions in bulk RNA-seq data.
![pipeline](./Figures/Workflow_new.png)

How to cite `locCSN`
-------------------
This work has not published yet, please see [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.02.13.431104v1).

Get Started
-----------------
First install `locCSN` pacakge. All locCSN python functions are also stored in [Python Folder](https://github.com/xuranw/locCSN/tree/main/Python). 


```
pip install locCSN
```

Please download datasets stored in [DataStore](https://github.com/xuranw/locCSN/tree/main/DataStore). 

## Dataset Summary

Dataset            |  Chutype          | ASD Brain              |  Brain Cortex Atlas
-------------------|-------------------|------------------------|---------------------------
Reference          |  Chu et al.(2016) | Velmeshev et al.(2019) | Polioudakis et al.(2019) 
\# of cell         |  1018             |  104,559               |    35,543                 
genes for analysis |  51               |     942                |     444                      
Data Availability  |    GSE75748       |  PRJNA434002           |  [website](http://solo.bmap.ucla.edu/shiny/webapp) 

## Simple Example of locCSN using Chutype Dataset

In this example, we reproduce the results of Chutype dataset in paper. 

### Load Datasets

There are 51 marker genes and 1018 cells from 7 cell types. The gene expression are stored in [logChumaker.txt](https://github.com/xuranw/locCSN/blob/main/DataStore/Chutype/logChumarker.txt) and corresponding cell types in [chutypectname.txt](https://github.com/xuranw/locCSN/blob/main/DataStore/Chutype/chutypectname.txt). Cell types are H1, H9, DEC, EC, HFF, NPC and TF. In our paper, we focus on cell type DEC and NPC.

```
# Import packages
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Set path to data
os.chdir('yourpathtodata/')

# read in Chutype dataset
data = sc.read_text('logChumarker.txt')
data.shape
data = data.transpose() # 1018 cells * 51 genes

cell_type = pd.read_csv('chutypectname.txt', sep = ' ')
data.obs = cell_type   # The observation are labeled with cell types.

# Plot the Heatmap of gene expression
sc.pl.heatmap(data, data.var.index, groupby= "cell_type", dendrogram = False, swap_axes = True, show_gene_labels= True, cmap='Wistia', figsize=(8,6))
```

![heatmap](./Figures/heatmapheat_exprs.png)


### Calculate Pearson's Correlation

After loading gene expression matrix and cell types, we first show the absolute Pearson's correlation for DEC and NPC cells. 
```
data_dec = data[data.obs.cell_type == "DEC", ]
X_dec = data_dec.X.transpose()
data_npc = data[data.obs.cell_type == 'NPC', ]
X_npc = data_npc.X.transpose()

corr_dec = np.corrcoef(X_dec)
corr_npc = np.corrcoef(X_npc)


np.fill_diagonal(corr_dec, 0)
np.fill_diagonal(corr_npc, 0)


plt.subplot(1, 2, 1)
plt.imshow(abs(corr_dec), vmin=0, vmax=0.7, cmap='RdPu')
plt.title('DEC', fontweight ="bold")
plt.subplot(1, 2, 2)
plt.imshow(abs(corr_npc), vmin=0, vmax=0.7, cmap='RdPu')
plt.title("NPC", fontweight = "bold")
plt.suptitle("Absolute Pearson`s Correlation", fontsize = 14, fontweight = "bold")
```
The heatmaps for absolute Pearson`s correlations is 

![abs_corr](./Figures/Abs_cor_2ct.png)


### Calculate CSN test statistics

Now we calculate the CSN test statistics using function `csn` for cell type DEC and NPC. 

```
import time
start = time.time()
csn_dec = locCSN.csn(X_dec, dev = True)
end = time.time()
print(end - start) 
start = time.time()
csn_npc = locCSN.csn(X_npc, dev = True)
end = time.time()
print(end - start) 
#1275 pairs need calculation
#60.3697772026062
#903 pairs need calculation
#35.72847938537598
```
Now we show what function `csn` produces.
```
type(csn_dec) 
# list
len(csn_dec) # 138 cells
# Let's see the test statistics for the first cell in DEC
plt.imshow(csn_dec[0].toarray(), vmin = -6, vmax = 6, cmap = 'coolwarm')
plt.title('DEC one cell', fontweight = "bold")
plt.colorbar()
#plt.savefig('dec_one_cell.png')

```
![one_cell](./Figures/dec_one_cell.png)

We compute each pair of genes and store test statistics in an upper diagnol matrix.


```
from scipy.stats import norm

# Cutoff at norm(0.95)
csn_mat = [(item > norm.ppf(0.95)).astype(int) for item in csn_dec]
avgcsn_dec = sum(csn_mat).toarray()/len(csn_mat) + np.transpose(sum(csn_mat).toarray()/len(csn_mat))
csn_mat = [(item > norm.ppf(0.95)).astype(int) for item in csn_npc]
avgcsn_npc = sum(csn_mat).toarray()/len(csn_mat) + np.transpose(sum(csn_mat).toarray()/len(csn_mat))

plt.subplot(1, 2, 1)
plt.imshow(avgcsn_dec, cmap = "Greens", vmin = 0, vmax = 0.7)
plt.title('DEC', fontweight ="bold")
plt.subplot(1, 2, 2)
plt.imshow(avgcsn_npc, cmap = "Greens", vmin = 0, vmax = 0.7)
plt.title('NPC', fontweight = 'bold')
plt.suptitle("Averaged CSN, cut at alpha = 0.05", fontsize=14, fontweight = "bold")
```
The heatmaps for DEC and NPC are 

![avgcsn](./Figures/Avg_csn_2ct.png)

## Comparing networks between two groups of cells

For comparison between two groups of cells using CSNs, we use the dataset from Brain Cortex Atlas. We focus on 942 expressed SFARI ASD genes. The comparison sLED, DISTp. sLED in this [GitHub repo](https://github.com/lingxuez/sLED)

### CSN 
I will provide the code for generating the CSN matrix, but the CSN results can be found in 

#### sLED comparison 

#### DISTp


## Trajectory analysis using Brain Cortex Atlas Dataset

PisCES is a Matlab package in this [GitHub repo](https://github.com/letitiaLiu/PisCES)



## References
* Dai, Hao, et al. "Cell-specific network constructed by single-cell RNA sequencing data." Nucleic acids research 47.11 (2019).
* Chu, Li-Fang, et al. "Single-cell RNA-seq reveals novel regulators of human embryonic stem cell differentiation to definitive endoderm." Genome biology 17.1 (2016).
* Velmeshev, Dmitry, et al. "Single-cell genomics identifies cell type-specific molecular changes in autism." Science 364.6441 (2019).
* Polioudakis, Damon, et al. "A single cell transcriptomic atlas of human neocortical development during mid-gestation." Neuron 103.5 (2019).
* Zhu, Lingxue, et al. "Testing high-dimensional covariance matrices, with application to detecting schizophrenia risk genes." The annals of applied statistics 11.3 (2017).
* Liu, Fuchen, et al. "Global spectral clustering in dynamic networks." Proceedings of the National Academy of Sciences 115.5 (2018).
