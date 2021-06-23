ESCO Simulations in Supplementary Note 2
================

Before execution, we need to load
[ESCO](https://github.com/JINJINT/ESCO) package.

``` r
library("devtools")
devtools::install_github("JINJINT/ESCO")

library(ESCO)
```

## Simulation for Compare CSN with BigSCale correlation

``` r
sim.dis = escoSimulateGroups(nGenes = 10000, nCells = 5000, 
                              group.prob = c(0.5, 0.5), lib.loc = 8, withcorr = TRUE,
                              dropout.mid = 3, alpha_mean = 0.05, depth_mean = 800, depth_sd = 10, 
                              de.facLoc = c(2, 1), corr = corr.true, trials = 1, verbose = TRUE)
#saveRDS(sim.dis, file = paste0(dir.data, 'sim_discrete_forbs_2.rds'))
genegroup = paste0("Group", rowData(sim.dis)$GeneGroup)
genegroup[which(genegroup=="Group0")] = "None"
geneinfo = data.frame(genes = rowData(sim.dis)$Gene, 
                      newcelltype = as.factor(genegroup))

# organize the cell info
cellinfo = data.frame(cells = colData(sim.dis)$Cell, 
                      newcelltype= as.factor(colData(sim.dis)$Group))

# get the data
datalist = list("simulated truth" = assays(sim.dis)$TrueCounts, 
                "zero-inflated" = assays(sim.dis)$counts, 
                "down-sampled" = assays(sim.dis)$observedcounts)
degeneinfo = geneinfo[which(geneinfo$newcelltype!="None"),]
degeneinfo$newcelltype = droplevels(degeneinfo$newcelltype)

heatdata(datalist, cellinfo = cellinfo, geneinfo = geneinfo, 
         size = 1, ncol = 3)
# plot GCN for all marker genes (i.e. DE genes) across all cell groups
degeneinfo = geneinfo[which(geneinfo$newcelltype!="None"),]
degeneinfo$newcelltype = droplevels(degeneinfo$newcelltype)
degcnlist = lapply(datalist, function(data)gcn(data, genes = degeneinfo$genes))
heatgcn(degcnlist, geneinfo = degeneinfo, size = 2, ncol = 3)

log.cpm = log(apply(datalist$`simulated truth`, 2, function(x){x/sum(x)})*10^6+1)
counts.marker = datalist$`simulated truth`[degeneinfo$genes, ]
log.cpm.marker = log.cpm[degeneinfo$genes, ]

log.marker.group1 = log.cpm.marker[, cellinfo$newcelltype == 'Group1']
log.marker.group2 = log.cpm.marker[, cellinfo$newcelltype == 'Group2']

cor.marker.group1 = cor(t(log.marker.group1))
cor.marker.group2 = cor(t(log.marker.group2))


cor.heat.disc = plot.cor.heat(log.cpm.marker, cellinfo$newcelltype, na2zer = T, axis.text = F, grid = F)

save(datalist, degeneinfo, cellinfo, file = paste0(dir.data, 'sim_dis_for_BigSCale.RData'))
```
