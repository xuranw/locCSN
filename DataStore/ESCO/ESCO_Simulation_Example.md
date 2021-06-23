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

## Simulation for true counts without technical noise
```r
corr.single = list()
corr.single[[1]] = diag(75)
#diag(corr.single[[1]]) = 1
#corr.single[[1]]
corr.single[[1]][1:15, 1:15] <- corr.single[[1]][16:30, 16:30] <- corr.single[[1]][31:45, 31:45] <- 
  corr.single[[1]][46:60, 46:60] <- corr.single[[1]][61:75, 61:75]  <- 0.9
corr.single[[1]][1:15, 16:30] <- corr.single[[1]][16:30, 31:45] <- corr.single[[1]][31:45, 46:60] <-
  corr.single[[1]][46:60, 61:75] <- corr.single[[1]][16:30, 1:15] <- corr.single[[1]][31:45, 16:30] <- 
  corr.single[[1]][46:60, 31:45] <- corr.single[[1]][61:75, 46:60] <- 0.7
corr.single[[1]][1:15, 31:45] <- corr.single[[1]][16:30, 46:60] <- corr.single[[1]][31:45, 61:75] <-
  corr.single[[1]][31:45, 1:15] <- corr.single[[1]][46:60, 16:30] <- 
  corr.single[[1]][61:75, 31:45] <- 0.5
corr.single[[1]][1:15, 46:60] <- corr.single[[1]][16:30, 61:75] <- corr.single[[1]][46:60, 1:15] <-
  corr.single[[1]][61:75, 16:30] <- 0.3
corr.single[[1]][1:15, 61:75] <- corr.single[[1]][61:75, 1:15] <- 0.1
diag(corr.single[[1]]) = 1

alpha = 0.1; lib.loc = 9;
sim.single <- escoSimulateSingle(nGenes = 100, nCells = 200,  withcorr = TRUE, corr = corr.single, verbose = TRUE)
corr.single = metadata(sim.single)$Params@corr

gene.order = colnames(corr.single[[1]])[hclust(dist(corr.single[[1]]))$order]
corr.single[[1]][gene.order[seq(1, 75, 15)], gene.order[seq(1, 75, 15)]]
gene.order = gene.order[c( 1:15, 16:30, 31:45, 61:75, 46:60)]
image(corr.single[[1]][gene.order, gene.order])

gene.order = c(gene.order, setdiff(paste0('Gene', 1:100), gene.order))

# get the data
datalist = list("simulated truth" = assays(sim.single)$TrueCounts, 
                "zero-inflated" = assays(sim.single)$counts, 
                "down-sampled" = assays(sim.single)$observedcounts)


#log.cpm = log(apply(datalist$`down-sampled`, 2, function(x){x/sum(x)})*10^6+1)
log.cpm = log(apply(datalist$`simulated truth`, 2, function(x){x/sum(x)})*10^6 + 1)
log.cpm = log.cpm[gene.order, ]
write.table(log.cpm, file = paste0(dir.save, 'simsingle_logcpm_trail.txt'))
```

## Simulation for true counts with technical noise
### Strong Connection
```r
corr.single = list()
corr.single[[1]] = diag(100)
corr.single[[1]][1:25, 1:25] <- corr.single[[1]][26:50, 26:50] <- 
  corr.single[[1]][51:75, 51:75] <- corr.single[[1]][76:100, 76:100] <- 0.95
diag(corr.single[[1]]) = 1

sim.single <- escoSimulateSingle(nGenes = 100, nCells = 200, lib.loc = lib.loc, 
                                 withcorr = TRUE, corr = corr.single, verbose = TRUE,
                                 dropout.type = 'downsample', alpha_mean = alpha)

mean(assays(sim.single)$observedcounts == 0)

corr.single = metadata(sim.single)$Params@corr
gene.order = paste0('Gene', hclust(dist(corr.single[[1]]))$order)
image(corr.single[[1]][gene.order, gene.order])

# get the data
datalist = list("simulated truth" = assays(sim.single)$TrueCounts, 
                "zero-inflated" = assays(sim.single)$counts, 
                "down-sampled" = assays(sim.single)$observedcounts)

log.cpm = log2(apply(datalist$`down-sampled`, 2, function(x){x/sum(x)})*10^6 + 1)

log.cpm = log.cpm[gene.order, ]
write.table(log.cpm, file = paste0(dir.save, 'simsingle_logcpm_trail.txt'))
save(sim.single, file = paste0(dir.save, 'sim_single.RData'))
```

### With Weak Connections
```r
corr.single = list()
corr.single[[1]] = diag(100)
corr.single[[1]][1:25, 1:25] <- corr.single[[1]][26:50, 26:50] <- 
  corr.single[[1]][51:75, 51:75] <- corr.single[[1]][76:100, 76:100] <- 0.95
corr.single[[1]][1:25, 26:50] <- corr.single[[1]][26:50, 51:75] <- corr.single[[1]][26:50, 1:25] <- 
  corr.single[[1]][51:75, 26:50] <- 0.5
diag(corr.single[[1]]) = 1

sim.single <- escoSimulateSingle(nGenes = 100, nCells = 200, lib.loc = lib.loc, 
                                 withcorr = TRUE, corr = corr.single, verbose = TRUE,
                                 dropout.type = 'downsample', alpha_mean = alpha)

mean(assays(sim.single)$observedcounts == 0)

corr.single = metadata(sim.single)$Params@corr
gene.order = paste0('Gene', hclust(dist(corr.single[[1]]))$order)
image(corr.single[[1]][gene.order, gene.order])

# get the data
datalist = list("simulated truth" = assays(sim.single)$TrueCounts, 
                "zero-inflated" = assays(sim.single)$counts, 
                "down-sampled" = assays(sim.single)$observedcounts)

log.cpm = log2(apply(datalist$`down-sampled`, 2, function(x){x/sum(x)})*10^6 + 1)

log.cpm = log.cpm[gene.order, ]
write.table(log.cpm, file = paste0(dir.save, 'simsingle_logcpm_trail.txt'))
save(sim.single, file = paste0(dir.save, 'sim_single.RData'))
```
