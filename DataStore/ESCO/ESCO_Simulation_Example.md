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
log.cpm = read.table(paste0(dir.save, 'simsingle_logcpm_trail.txt'))

cor.pearson = cor(t(log.cpm)); diag(cor.pearson) = NA;
m.cor.temp = data.frame(cor = c(cor.pearson), X = factor(rep(gene.order, 100), levels = gene.order), 
                        Y = factor(rep(gene.order, each = 100), levels = gene.order))
m.cor.temp$cor[is.na(m.cor.temp$cor)] = 0

jpeg(paste0(dir.save, 'trial_single_pearson.jpeg'), width = 300, height = 300)
ggplot(m.cor.temp, aes(X, Y, fill = cor)) + geom_tile() + scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1)) + 
  theme_xuran() + ggtitle('Pearson\'s Correlation') + geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5) + coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9') + geom_text(x = 22.5, y = 7.5, label = '0.7') + 
  geom_text(x = 37.5, y = 7.5, label = '0.5') + geom_text(x = 52.5, y = 7.5, label = '0.3') + 
  geom_text(x = 67.5, y = 7.5, label = '0.1') + geom_text(x = 82.5, y = 7.5, label = '0') + 
  theme(axis.text = element_blank(), axis.text.x = element_blank())
dev.off()

cor.true = diag(100); cor.true[1:75, 1:75] = corr.single[[1]][gene.order[1:75], gene.order[1:75]]
m.cor.true = data.frame(rho = c(cor.true), X = factor(rep(gene.order, 100), levels = gene.order), 
                        Y = factor(rep(gene.order, each = 100), levels = gene.order))

#jpeg(paste0(dir.save, 'trial_single_true.jpeg'), width = 300, height = 300)
pdf(paste0(dir.save, 'trial_single_true.pdf'), width = 6, height = 6)
p_true = ggplot(m.cor.true, aes(X, Y, fill = rho)) + geom_tile() + 
  scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1), name = expression(rho)) + 
  theme_minimal() + ggtitle('True Connection') + 
  geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9', size = 3) + 
  geom_text(x = 22.5, y = 7.5, label = '0.7', size = 3) + 
  geom_text(x = 37.5, y = 7.5, label = '0.5', size = 3) + 
  geom_text(x = 52.5, y = 7.5, label = '0.3', size = 3) + 
  geom_text(x = 67.5, y = 7.5, label = '0.1', size = 3) + 
  geom_text(x = 82.5, y = 7.5, label = '0', size = 3) + 
  theme(axis.text = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5))
dev.off()

rho.mat

rho.flat = NULL
for(i in 1:99){
  for(j in (i+1):100){
    if(i <= 75 & j <= 75){
      rho.flat = c(rho.flat, max(0.9 - abs(floor((i-1)/15) - floor((j-1)/15))*0.2, 0))
    }else{
      rho.flat = c(rho.flat, 0)
    }
  }
}
rho.flat = round(rho.flat, 1)

rho.mat = diag(100)
k = 0
for(i in 1:99){
  for(j in (i+1):100){
    k = k+1
    rho.mat[i, j] <- rho.mat[j, i] <- rho.flat[k]
  }
}
############# Quantile window sizes ################
csn.o005.flat = data.matrix(read.table(paste0(dir.save, 'csn_o005_flat.txt')))
csn.o01.flat = data.matrix(read.table(paste0(dir.save, 'csn_o01_flat.txt')))
csn.o015.flat = data.matrix(read.table(paste0(dir.save, 'csn_o015_flat.txt')))
csn.o02.flat = data.matrix(read.table(paste0(dir.save, 'csn_o02_flat.txt')))

rho = 0
csn.rho.0 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                  c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))
rho = 0.1
csn.rho.01 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                   c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))

rho = 0.3
csn.rho.03 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                   c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))

rho = 0.5
csn.rho.05 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                   c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))

rho = 0.7
csn.rho.07 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                   c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))

rho = 0.9
csn.rho.09 = cbind(c(csn.o005.flat[rho.flat == rho, ]), c(csn.o01.flat[rho.flat == rho, ]), 
                   c(csn.o015.flat[rho.flat == rho, ]), c(csn.o02.flat[rho.flat == rho, ]))


table(rowSums(csn.rho.0 == 0))
#0      1      2      3      4 
#255577   1487   1480   2714 395742 
# this means part of the gene 
csn.rho.small.0 = csn.rho.0[rowSums(csn.rho.0 == 0) < 4, ]
csn.rho.small.01 = csn.rho.01[rowSums(csn.rho.01 == 0) < 4, ]
csn.rho.small.03 = csn.rho.03[rowSums(csn.rho.03 == 0) < 4, ]
csn.rho.small.05 = csn.rho.05[rowSums(csn.rho.05 == 0) < 4, ]
csn.rho.small.07 = csn.rho.07[rowSums(csn.rho.07 == 0) < 4, ]
csn.rho.small.09 = csn.rho.09[rowSums(csn.rho.09 == 0) < 4, ]

q.name = c('quantile 5%', 'quantile 10%', 'quantile 15%', 'quantile 20%')

list.csn.rho.small = list(csn.rho.small.0, csn.rho.small.01, csn.rho.small.03, 
                          csn.rho.small.05, csn.rho.small.07, csn.rho.small.09)
rho.name = paste0('rho = ', c(0, seq(0.1, 0.9, 0.2)))

library(plyr)

my_label_parsed <- function (variable, value) {
  if (variable == "quantile") {
    return(as.character(value))
  } else {
    llply(as.character(value), function(x) parse(text = x))    
  }
}

rho.label = c("rho == 0", "rho == 0.1", "rho == 0.3", "rho == 0.5", "rho == 0.7", "rho == 0.9")
m.csn.rho.small = NULL
for(i in 1:length(rho.name)){
  temp = data.frame(test.stat = c(list.csn.rho.small[[i]]), 
                    quantile = factor(rep(q.name, each = nrow(list.csn.rho.small[[i]])), levels = q.name), 
                    rho = rep(rho.name[i], nrow(list.csn.rho.small[[i]])*4), 
                    rho.label = rep(rho.label[i], nrow(list.csn.rho.small[[i]])*4))
  
  m.csn.rho.small = rbind(m.csn.rho.small, temp)
}

#jpeg(filename = paste0(dir.save, 'hist_qtl_trail.jpg'), width = 400, height = 300)
pdf(paste0(dir.save, 'hist_qtl_trail.pdf'), width = 7.5, height = 6)
p_hist_q = ggplot(m.csn.rho.small, aes(test.stat)) + geom_histogram(bins = 20) + 
  facet_grid(rho.label~quantile, scales = 'free_y', switch = "y", labeller = my_label_parsed) + 
  theme_minimal() + ggtitle('Quantile Window Sizes') + 
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.3), strip.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8))
dev.off()
t(apply(csn.rho.small.0, 2, function(x){quantile(x, c(0.9, 0.95, 0.99, 1))}))
# 90%      95%      99%     100%
#[1,] 3.250181 3.250181 4.842438 12.15366
#[2,] 2.097809 2.848367 4.912354 11.19207
#[3,] 2.259124 2.867098 5.490291 11.09749
#[4,] 2.419285 3.287027 6.173242 11.34821

t(apply(csn.rho.small.01, 2, function(x){quantile(x, c(0.9, 0.95, 0.99, 1))}))

## ROC curve: assign rho = 0, 0.1, and 0.3 with no connections for all cells 
## and rho > 0.3 with connection for all cells.
# 50% of all samples be TN, discard 0 < rho < 0.3, >= 0.3
sum(rho.flat > 0.3)  # 2100
selected.pairs = sort(c(sample(which(rho.flat == 0), sum(rho.flat > 0.3)), which(rho.flat > 0.3)))

bench.rho = rho.flat
bench.rho[-selected.pairs] = NA
bench.rho[which(bench.rho > 0)] = 1
bench.mat = diag(100)
k = 0
for(i in 1:99){
  for(j in (i+1):100){
    k = k+1
    bench.mat[i, j] <- bench.mat[j, i] <- bench.rho[k]
  }
}
m.bench.temp = data.frame(benchmark = c(bench.mat), X = factor(rep(gene.order, 100), levels = gene.order), 
                        Y = factor(rep(gene.order, each = 100), levels = gene.order))
#jpeg(paste0(dir.save, 'trial_single_bench.jpg'), width = 400, height = 400)
pdf(paste0(dir.save, 'trial_single_bench.pdf'), width = 8, height = 8)
p_bench = ggplot(m.bench.temp, aes(X, Y, fill = benchmark)) + geom_tile() + scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1)) + 
  theme_minimal() + ggtitle('Benchmark Matrix') + 
  geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9', size = 3) + 
  geom_text(x = 22.5, y = 7.5, label = '0.7', size = 3) + 
  geom_text(x = 37.5, y = 7.5, label = '0.5', size = 3) + 
  geom_text(x = 52.5, y = 7.5, label = '0.3', size = 3) + 
  geom_text(x = 67.5, y = 7.5, label = '0.1', size = 3) + 
  geom_text(x = 82.5, y = 7.5, label = '0', size = 3) +   
  theme(axis.text = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5))
dev.off()

connect.rho.temp = as.matrix(rho.flat[selected.pairs] > 0)

list.csn.rho.select = list(csn.o005.flat[selected.pairs, ], csn.o01.flat[selected.pairs, ], 
                    csn.o015.flat[selected.pairs, ], csn.o02.flat[selected.pairs, ])
connect.rho = connect.rho.temp[, rep(1, 200)]

q.name = c('quantile 5%', 'quantile 10%', 'quantile 15%', 'quantile 20%')
lapply(list.csn.rho.select, range)

t.list = c(seq(-10, -0.5, 0.5), seq(0, 12.6, 0.05))
T.pos = NULL; F.pos = NULL; Acc = NULL;
for(i in 1:4){
  t.pos = NULL
  f.pos = NULL
  acc = NULL
  for(t in t.list){
    t.pos = c(t.pos, sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 11)/sum(connect.rho == 1) )  # True positive
    f.pos = c(f.pos, sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 1)/ sum(connect.rho == 0) )  # False positive
    acc = c(acc, (sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 11)+sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 0))/(nrow(connect.rho)*ncol(connect.rho)) )
  }
  T.pos = cbind(T.pos, t.pos)
  F.pos = cbind(F.pos, f.pos)
  Acc = cbind(Acc, acc)
}
m.roc = data.frame(True.positive = c(T.pos), False.positive = c(F.pos),# threshold = 
                     window = factor(rep(q.name, each = length(t.list)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ROC_qtl.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ROC_qtl.pdf'), width = 7, height = 5)
p_roc_q = ggplot(m.roc, aes(False.positive, True.positive, col = window)) + geom_line(size = 0.5) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + scale_color_brewer(palette = 'Set1')+ 
  theme_minimal() + theme(legend.position = 'none', axis.text = element_text(size = 8)) + ggtitle('ROC')
dev.off()

jpeg(filename = paste0(dir.save, 'ROC_qtl_s.jpg'), width = 800, height = 250)
pdf(paste0(dir.save, 'ROC_qtl_s.pdf'), width = 12, height = 8)
ggplot(m.roc, aes(False.positive, True.positive, col = window)) + geom_line(size = 0.5) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + scale_color_brewer(palette = 'Set1')+ theme_xuran() + 
  facet_wrap( ~ window, nrow = 1) + ggtitle('ROC')
dev.off()

m.acc = data.frame(ACC = c(Acc), threshold = rep(t.list, 4), 
                   window = factor(rep(q.name, each = length(t.list)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ACC_qtl_s.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ACC_qtl_s.pdf'), width = 6, height = 4)
p_acc_q = ggplot(m.acc, aes(threshold, ACC, col = window)) + geom_line(size = 0.5) + #coord_fixed(ratio = 250) +
  scale_color_brewer(palette = 'Set1') + ggtitle('ACC') + theme_minimal() + 
  theme(legend.position = 'none', axis.text = element_text(size = 8), aspect.ratio=1/1)
dev.off()

Acc.q = NULL
t.q.list = qnorm(seq(0.85, 1-0.0001, 0.0001))
for(i in 1:4){
  acc = NULL
  for(t in t.q.list){
    acc = c(acc, (sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 11)+
                    sum((list.csn.rho.select[[i]] > t) + 10*connect.rho == 0))/(nrow(connect.rho)*ncol(connect.rho)) )
  }
  Acc.q = cbind(Acc.q, acc)
}


#max accuracy is around 0.9
m.acc.q = data.frame(ACC = c(Acc.q), threshold = rep(seq(0.85, 1-0.0001, 0.0001), 4), 
                   window = factor(rep(q.name, each = length(t.q.list)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ACC_qtl_q.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ACC_qtl_q.pdf'), width = 6, height = 4)
p_acc_q_q = ggplot(m.acc.q, aes(threshold, ACC, col = window)) + geom_line(size = 0.5) + 
  #coord_fixed(ratio = 2.5) + 
  scale_color_brewer(palette = 'Set1') + ggtitle('ACC') + 
  theme_minimal() + xlab(expression(paste("threshold (1-", alpha,")"))) + 
  theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8), aspect.ratio=1/1)
dev.off()


index_mat = NULL
for(i in 1:(G-1)){
  for(j in (i+1):G){
    index_mat = rbind(index_mat, c(i, j))
  }
}
list.csn.rho = list(csn.o005.flat, csn.o01.flat, csn.o015.flat, csn.o02.flat)

 #t.s = t.list[182]
G = length(gene.order); N = ncol(csnflat)
m.avgcsn = NULL
for(i in 1:4){
  csnflat = list.csn.rho[[i]];
  avgcsn.temp = data.frame(avgcsn = rep(rowMeans(csnflat > qnorm(0.95)), 2), 
                           X = factor(gene.order[c(index_mat[, 1], index_mat[, 2])], levels = gene.order), 
                           Y = factor(gene.order[c(index_mat[, 2], index_mat[, 1])], levels = gene.order))
  avgcsn.temp = rbind(avgcsn.temp, data.frame(avgcsn = rep(1, G), X = gene.order, Y = gene.order))
  avgcsn.temp$window = rep(q.name[i],  nrow(avgcsn.temp))                 
  m.avgcsn = rbind(m.avgcsn, avgcsn.temp)
}
m.avgcsn$window = factor(m.avgcsn$window, levels = q.name)

#jpeg(filename = paste0(dir.save, 'heatmap_qtl_qnorm.jpg'), width = 800, height = 250)
pdf(paste0(dir.save, 'heatmap_qtl_qnorm.pdf'), width = 12, height = 5)
p_avgcsn_q = ggplot(m.avgcsn, aes(X, Y, fill = avgcsn)) + geom_tile() + scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1)) + 
  theme_minimal() + facet_wrap(~ window, nrow = 1) + 
  ggtitle(expression(paste('Averaged CSN from Quantile Windows, cut at ', alpha, '= 0.05'))) + 
  geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9', size = 3) + 
  geom_text(x = 22.5, y = 7.5, label = '0.7', size = 3) + 
  geom_text(x = 37.5, y = 7.5, label = '0.5', size = 3) + 
  geom_text(x = 52.5, y = 7.5, label = '0.3', size = 3) + 
  geom_text(x = 67.5, y = 7.5, label = '0.1', size = 3) + 
  geom_text(x = 82.5, y = 7.5, label = '0', size = 3) + 
  theme(axis.text = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5))
dev.off()

G = length(gene.order); N = ncol(csnflat)
m.avgcsn.new = NULL
for(i in 1:4){
  csnflat = list.csn.rho[[i]];
  avgcsn.temp = data.frame(avgcsn = rep(rowMeans(csnflat > qnorm(0.99)), 2), 
                           X = factor(gene.order[c(index_mat[, 1], index_mat[, 2])], levels = gene.order), 
                           Y = factor(gene.order[c(index_mat[, 2], index_mat[, 1])], levels = gene.order))
  avgcsn.temp = rbind(avgcsn.temp, data.frame(avgcsn = rep(NA, G), X = gene.order, Y = gene.order))
  avgcsn.temp$window = rep(q.name[i],  nrow(avgcsn.temp))                 
  m.avgcsn.new = rbind(m.avgcsn.new, avgcsn.temp)
}
m.avgcsn.new$window = factor(m.avgcsn.new$window, levels = q.name)

jpeg(filename = paste0(dir.save, 'heatmap_qtl_ts.jpg'), width = 800, height = 250)
p_avgcsn_q_new = ggplot(m.avgcsn.new, aes(X, Y, fill = avgcsn)) + geom_tile() + scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1)) + 
  theme_minimal() + facet_wrap(~ window, nrow = 1) + ggtitle(expression(paste('Averaged CSN, cut at ', alpha, '= 0.01'))) + 
  geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9', size = 3) + 
  geom_text(x = 22.5, y = 7.5, label = '0.7', size = 3) + 
  geom_text(x = 37.5, y = 7.5, label = '0.5', size = 3) + 
  geom_text(x = 52.5, y = 7.5, label = '0.3', size = 3) + 
  geom_text(x = 67.5, y = 7.5, label = '0.1', size = 3) + 
  geom_text(x = 82.5, y = 7.5, label = '0', size = 3) + 
  theme(axis.text = element_blank(), axis.text.x = element_blank())
dev.off()



################ Standard deviation window size ############
csn.d005.flat = data.matrix(read.table(paste0(dir.save, 'csn_d005_flat.txt')))
csn.d01.flat = data.matrix(read.table(paste0(dir.save, 'csn_d01_flat.txt')))
csn.d015.flat = data.matrix(read.table(paste0(dir.save, 'csn_d015_flat.txt')))
csn.d02.flat = data.matrix(read.table(paste0(dir.save, 'csn_d02_flat.txt')))

rho = 0
csn.rho.d0 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                  c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))
rho = 0.1
csn.rho.d01 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                   c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))

rho = 0.3
csn.rho.d03 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                   c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))

rho = 0.5
csn.rho.d05 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                   c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))

rho = 0.7
csn.rho.d07 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                   c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))

rho = 0.9
csn.rho.d09 = cbind(c(csn.d005.flat[rho.flat == rho, ]), c(csn.d01.flat[rho.flat == rho, ]), 
                   c(csn.d015.flat[rho.flat == rho, ]), c(csn.d02.flat[rho.flat == rho, ]))


table(rowSums(csn.rho.d0 == 0))
#0      1      2      3      4 
#260604   1689     98     33 394576 

# this means part of the gene 
csn.rho.small.d0 = csn.rho.d0[rowSums(csn.rho.d0 == 0) < 4, ]
csn.rho.small.d01 = csn.rho.d01[rowSums(csn.rho.d01 == 0) < 4, ]
csn.rho.small.d03 = csn.rho.d03[rowSums(csn.rho.d03 == 0) < 4, ]
csn.rho.small.d05 = csn.rho.d05[rowSums(csn.rho.d05 == 0) < 4, ]
csn.rho.small.d07 = csn.rho.d07[rowSums(csn.rho.d07 == 0) < 4, ]
csn.rho.small.d09 = csn.rho.d09[rowSums(csn.rho.d09 == 0) < 4, ]
q.name = c('SD 5%', 'SD 10%', 'SD 15%', 'SD 20%')
list.csn.rho.small = list(csn.rho.small.d0, csn.rho.small.d01, csn.rho.small.d03, 
                          csn.rho.small.d05, csn.rho.small.d07, csn.rho.small.d09)
rho.name = paste0('rho = ', c(0, seq(0.1, 0.9, 0.2)))
m.csn.rho.small.d = NULL
for(i in 1:length(rho.name)){
  temp = data.frame(test.stat = c(list.csn.rho.small[[i]]), 
                    quantile = factor(rep(q.name, each = nrow(list.csn.rho.small[[i]])), levels = q.name), 
                    rho = rep(rho.name[i], nrow(list.csn.rho.small[[i]])*4), 
                    rho.label = rep(rho.label[i], nrow(list.csn.rho.small[[i]])*4))
  
  m.csn.rho.small.d = rbind(m.csn.rho.small.d, temp)
}
jpeg(filename = paste0(dir.save, 'hist_sd_trail.jpg'), width = 400, height = 300)
#pdf(paste0(dir.save, 'hist_sd_trail.pdf'), width = 7.5, height = 6)
p_hist_sd = ggplot(m.csn.rho.small.d, aes(test.stat)) + geom_histogram(bins = 20) + 
  facet_grid(rho.label ~ quantile, scales = 'free_y', switch = "y", labeller = my_label_parsed) + 
  theme_minimal() + ggtitle('Standard Deviation Window Sizes') + 
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.3), strip.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8))
dev.off()

list.csn.rho.d = list(csn.d005.flat, csn.d01.flat, csn.d015.flat, csn.d02.flat)
q.name = c('SD 5%', 'SD 10%', 'SD 15%', 'SD 20%')

list.csn.rho.d.select = lapply(list.csn.rho.d, function(x){x[selected.pairs, ]})
lapply(list.csn.rho.d.select, range)
t.list.d = c(seq(-11.5, -0.5, 0.5), seq(0, 14.2, 0.05))
T.pos.d = NULL; F.pos.d = NULL; Acc.d = NULL;
for(i in 1:4){
  t.pos = NULL
  f.pos = NULL
  acc = NULL
  for(t in t.list.d){
    t.pos = c(t.pos, sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 11)/sum(connect.rho == 1) )  # True positive
    f.pos = c(f.pos, sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 1)/ sum(connect.rho == 0) )  # False positive
    acc = c(acc, (sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 11)+
                    sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 0))/(nrow(connect.rho)*ncol(connect.rho)) )
  }
  T.pos.d = cbind(T.pos.d, t.pos)
  F.pos.d = cbind(F.pos.d, f.pos)
  Acc.d = cbind(Acc.d, acc)
}

m.roc.d = data.frame(True.positive = c(T.pos.d), False.positive = c(F.pos.d), 
                     window = factor(rep(q.name, each = length(t.list.d)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ROC_sd.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ROC_sd.pdf'), width = 7, height = 5)
p_roc_sd = ggplot(m.roc.d, aes(False.positive, True.positive, col = window)) + geom_line(size = 0.5) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + scale_color_brewer(palette = 'Set1')+ theme_minimal() + 
  theme(legend.position = 'none', axis.text = element_text(size = 8)) +
  ggtitle('ROC')
dev.off()

#jpeg(filename = paste0(dir.save, 'ROC_sd_s.jpg'), width = 800, height = 250)
pdf(paste0(dir.save, 'ROC_sd_s.pdf'), width = 12, height = 5)
ggplot(m.roc.d, aes(False.positive, True.positive, col = window)) + geom_line(size = 0.5) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + scale_color_brewer(palette = 'Set1')+ theme_xuran() + 
  facet_wrap( ~ window, nrow = 1) + ggtitle('ROC')
dev.off()

m.acc.d = data.frame(ACC = c(Acc.d), threshold = rep(t.list.d, 4), 
                   window = factor(rep(q.name, each = length(t.list.d)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ACC_sd_s.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ACC_sd_s.pdf'), width = 6, height = 4)
p_acc_sd = ggplot(m.acc.d, aes(threshold, ACC, col = window)) + geom_line(size = 0.5) + 
  #coord_fixed(ratio = 220) + 
  scale_color_brewer(palette = 'Set1') + ggtitle('ACC') + 
  theme_minimal() + theme(legend.position = 'none', axis.text = element_text(size = 8), aspect.ratio=1/1)
dev.off()

Acc.d.q = NULL
t.q.list = qnorm(seq(0.8, 1-0.0001, 0.0001))
for(i in 1:4){
  acc = NULL
  for(t in t.q.list){
    acc = c(acc, (sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 11)+
                    sum((list.csn.rho.d.select[[i]] > t) + 10*connect.rho == 0))/(nrow(connect.rho)*ncol(connect.rho)) )
  }
  Acc.d.q = cbind(Acc.d.q, acc)
}
m.acc.d.q = data.frame(ACC = c(Acc.d.q), threshold = rep(seq(0.8, 1-0.0001, 0.0001), 4), 
                     window = factor(rep(q.name, each = length(t.q.list)), levels = q.name))

#jpeg(filename = paste0(dir.save, 'ACC_sd_q.jpg'), width = 350, height = 250)
pdf(paste0(dir.save, 'ACC_sd_q.pdf'), width = 6, height = 4)
p_acc_sd_q = ggplot(m.acc.d.q, aes(threshold, ACC, col = window)) + geom_line(size = 0.5) + 
  #coord_fixed(ratio = 3.2) + 
  scale_color_brewer(palette = 'Set1') + ggtitle('ACC') + 
  theme_minimal() + xlab(expression(paste("threshold (1-", alpha,")"))) + 
  theme(legend.text = element_text(size = 8), axis.text = element_text(size = 8), aspect.ratio=1/1)
dev.off()


G = length(gene.order); N = ncol(csnflat)
m.avgcsn.d = NULL
for(i in 1:4){
  csnflat = list.csn.rho.d[[i]];
  avgcsn.temp = data.frame(avgcsn = rep(rowMeans(csnflat > qnorm(0.95)), 2), 
                           X = factor(gene.order[c(index_mat[, 1], index_mat[, 2])], levels = gene.order), 
                           Y = factor(gene.order[c(index_mat[, 2], index_mat[, 1])], levels = gene.order))
  avgcsn.temp = rbind(avgcsn.temp, data.frame(avgcsn = rep(1, G), X = gene.order, Y = gene.order))
  avgcsn.temp$window = rep(q.name[i],  nrow(avgcsn.temp))                 
  m.avgcsn.d = rbind(m.avgcsn.d, avgcsn.temp)
}
m.avgcsn.d$window = factor(m.avgcsn.d$window, levels = q.name)

#jpeg(filename = paste0(dir.save, 'heatmap_sd_qnorm.jpg'), width = 800, height = 250)
pdf(paste0(dir.save, 'heatmap_sd_qnorm.pdf'), width = 12, height = 5)
p_avgcsn_sd = ggplot(m.avgcsn.d, aes(X, Y, fill = avgcsn)) + geom_tile() + 
  scale_fill_distiller(palette = 'RdBu', limits = c(-1, 1)) + 
  theme_minimal() + facet_wrap(~ window, nrow = 1) + 
  ggtitle(expression(paste('Averaged CSN from SD Windows, cut at ', alpha, '= 0.05'))) + 
  geom_vline(xintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  geom_hline(yintercept = c(15, 30, 45, 60, 75)+0.5, size = 0.1) + 
  coord_fixed() + 
  geom_text(x = 7.5, y = 7.5, label = '0.9', size = 3) + 
  geom_text(x = 22.5, y = 7.5, label = '0.7', size = 3) + 
  geom_text(x = 37.5, y = 7.5, label = '0.5', size = 3) + 
  geom_text(x = 52.5, y = 7.5, label = '0.3', size = 3) + 
  geom_text(x = 67.5, y = 7.5, label = '0.1', size = 3) + 
  geom_text(x = 82.5, y = 7.5, label = '0', size = 3) + 
  theme(axis.text = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5))
dev.off()


m.roc.d$md = rep('sd', nrow(m.roc.d))
m.roc$md = rep('qtl', nrow(m.roc))
m.roc.full = rbind(m.roc.d, m.roc)
jpeg(filename = paste0(dir.save, 'ROC_sd_qtl_s.jpg'), width = 800, height = 500)
ggplot(m.roc.full, aes(False.positive, True.positive, col = window)) + geom_line(size = 0.5) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + scale_color_brewer(palette = 'Set1') + facet_wrap(~ window, nrow = 2) + 
  theme_xuran() + ggtitle('ROC')
dev.off()

plot1 <- plot_grid(p_true, p_bench, p_hist_q, p_hist_sd, align = "hv", nrow = 2, ncol = 2, 
                     labels = c("a", "b", "c", "d"), axis = "bt", rel_heights = c(1, 1.3))

ggsave(filename = paste0(dir.save, "simtruth_plot1.pdf"), plot = plot1, device = "pdf", 
       scale = 1, width = 8, height = 6)
ggsave(filename = paste0(dir.save, "simtruth_plot1.png"), plot = plot1, device = "png", 
       scale = 1, width = 8, height = 6)

top_row <- plot_grid(p_roc_q, p_acc_q, p_acc_q_q, p_roc_sd, p_acc_sd, p_acc_sd_q, align = 'hv', 
                        nrow = 2, ncol = 3, labels = c("a", "b", "c", "d", "e", "f"), 
                        rel_widths = c(1, 1, 1.5), axis = "bt")
bottom_row <-plot_grid(p_avgcsn_q, p_avgcsn_sd, align = 'hv', nrow = 2, ncol = 1, labels = c("g", "h"), 
                       axis = "bt")
plot2 <- plot_grid(top_row, bottom_row, nrow = 2, ncol = 1, labels = c("", ""), axis = "bt")
ggsave(filename = paste0(dir.save, "simtruth_plot2.pdf"), plot = plot2, device = "pdf", 
       scale = 1, width = 8, height = 10)
ggsave(filename = paste0(dir.save, "simtruth_plot2.png"), plot = plot2, device = "png", 
       scale = 1, width = 8, height = 10)


final_plot = plot_grid(plot1, plot2, align = "hv", nrow = 2, ncol = 1, 
                       labels = c("", ""), rel_heights =c(2.3, 4))
ggsave(filename = paste0(dir.save, "simtruth_combo.pdf"), plot = final_plot, device = "pdf", 
       scale = 1, width = 8, height = 15)
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
