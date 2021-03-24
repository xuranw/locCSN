
### Modify sLED for CSNs
library(sLED)

sLED.csn = function(X, Y, adj.beta = -1, rho = 1000, sumabs.seq = 0.2,
                    npermute = 100, useMC = FALSE, mc.cores = 1, seeds = NULL,
                    verbose = TRUE, niter = 20, trace = FALSE){
  D.hat = getDiffMatrix.csn(X, Y)
  pma.results <- sLEDTestStat(Dmat = D.hat, rho = rho, sumabs.seq = sumabs.seq,
                              niter = niter, trace = trace)
  Tn <- pma.results$stats
  n1 <- ncol(X)
  n2 <- ncol(Y)
  Z <- cbind(X, Y)
  permute.results <- sLEDpermute.csn(Z = Z, n1 = n1, n2 = n2, adj.beta = adj.beta,
                                     rho = rho, sumabs.seq = sumabs.seq, npermute = npermute,
                                     useMC = useMC, mc.cores = mc.cores, seeds = seeds, verbose = verbose,
                                     niter = niter, trace = trace)
  pVal = rowSums(permute.results$Tn.permute > Tn)/npermute
  return(c(pma.results, permute.results, list(Tn = Tn, pVal = pVal)))
}
getDiffMatrix.csn = function(X, Y){
  avg.x = rowMeans(X); avg.y = rowMeans(Y); avg.d = avg.x - avg.y;
  ng = ceiling(sqrt(length(avg.d)*2));
  D.hat = matrix(0, ng, ng); k = 0
  for(i in 1:(ng-1)){
    for(j in (i+1):ng){
      k = k+1
      D.hat[i, j] <- avg.d[k]
      D.hat[j, i] <- avg.d[k]
    }
  }
  return(D.hat)
}
sLEDTestStat <- function(Dmat, rho=1000, sumabs.seq=0.2,
                         niter=20, trace=FALSE) {
  ndim <- 1 ## only consider the first sparse eigenvector
  p <- ncol(Dmat)
  ntest <- length(sumabs.seq)

  results <- list()
  results$sumabs.seq <- sumabs.seq
  results$rho <- rho

  results$stats <- rep(NA, ntest)
  results$sign <- rep(NA, ntest)
  results$v <- matrix(NA, nrow=ntest, ncol=p)
  results$leverage <- matrix(NA, nrow=ntest, ncol=p)

  ## for each sparsity parameter
  for (i in 1:ntest) {
    sumabs <- sumabs.seq[i]
    pos.out <- symmPMD(Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)
    neg.out <- symmPMD(- Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)

    if (pos.out$d >= neg.out$d) {
      results$sign[i] <- "pos"
      results$stats[i] <- pos.out$d - rho
      results$v[i, ] <- pos.out$v
      results$leverage[i, ] <- (pos.out$v)^2
    } else {
      results$sign[i] <- "neg"
      results$stats[i] <- neg.out$d - rho
      results$v[i, ] <- neg.out$v
      results$leverage[i, ] <- (neg.out$v)^2
    }
  }

  return(results)
}
sLEDpermute.csn = function(Z, n1, n2, adj.beta=-1, rho=1000,
                           sumabs.seq=0.2, npermute=100,
                           useMC=FALSE, mc.cores=1, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE){
  if (verbose) {
    cat(npermute, "permutation started:\n")
  }

  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    perm.results <- parallel::mclapply(1:npermute, sLEDOnePermute.csn,
                                       Z=Z, n1=n1, n2=n2, seeds=seeds,
                                       sumabs.seq=sumabs.seq,
                                       adj.beta=adj.beta, rho=rho,
                                       verbose=FALSE, niter=niter, trace=trace,
                                       mc.cores=mc.cores, mc.preschedule=TRUE)
    } else {
      ## without parallelization
      perm.results <- lapply(1:npermute, sLEDOnePermute.csn,
                             Z=Z, n1=n1, n2=n2, seeds=seeds,
                             sumabs.seq=sumabs.seq,
                             adj.beta=adj.beta, rho=rho,
                             verbose=verbose, niter=niter, trace=trace)
    }

  ## extract test statistics and signs
  ntest <- length(sumabs.seq)
  Tn.permute <- matrix(NA, nrow=ntest, ncol=npermute)
  Tn.permute.sign <- matrix(NA, nrow=ntest, ncol=npermute)
  for (i in 1:npermute) {
    Tn.permute[, i] <- perm.results[[i]]$Tn.permute
    Tn.permute.sign[, i] <- perm.results[[i]]$Tn.permute.sign
  }

  if (verbose) {
    cat("permutations finished.", fill=TRUE)
  }

  return(list(Tn.permute = Tn.permute, Tn.permute.sign = Tn.permute.sign))
}
permuteIndex <- function(n1, n2){
  i.sample <- sample(1:(n1+n2), replace=FALSE)
  return(list(i1=i.sample[1:n1], i2=i.sample[(n1+1):(n1+n2)]))
}
sLEDOnePermute.csn = function(i, Z, n1, n2, seeds=NULL, sumabs.seq=0.2, adj.beta=-1, rho=1000,
                              verbose=TRUE, niter=20, trace=FALSE) {
  if (!is.null(seeds)){
    set.seed(seeds[i])
  }

  i.permute <- permuteIndex(n1, n2)
  D.permute <- getDiffMatrix.csn(X=Z[, i.permute$i1], Y=Z[, i.permute$i2])
  test.permute <- sLEDTestStat(Dmat=D.permute, rho=rho,
                               sumabs.seq=sumabs.seq, niter=niter, trace=trace)

  if (verbose  && (i %% 10)==0){
    cat(i, ",")
  }
  return(list(Tn.permute=test.permute$stats,
              Tn.permute.sign=test.permute$sign))
}

flat.to.matrix = function(flat.temp){
  ng = ceiling(sqrt(length(flat.temp)*2))
  D.mat = matrix(0, ng, ng); k = 0;
  for(i in 1:(ng-1)){
    for(j in (i+1):ng){
      k = k+1;
      D.mat[i, j] = flat.temp[k]
      D.mat[j, i] = flat.temp[k]
    }
  }
  return(D.mat)
}

flat.subset = function(remove.id, N){
  id.mat = matrix(0, N, N); k = 0;
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      k = k+1;
      id.mat[i, j] = k
    }
  }
  id.mat.sub = id.mat[-remove.id, -remove.id]
  ng = nrow(id.mat.sub)
  flat.id = rep(0, (ng-1)*(ng-2)/2); k = 0;
  for(i in 1:(ng-1)){
    for(j in (i+1):ng){
      k = k + 1;
      flat.id[k] = id.mat.sub[i, j]
    }
  }
  return(flat.id)
}

cor.nonzero = function(X){
  ng = ncol(X)
  n = nrow(X)
  A = diag(ng)
  for(i in 1:(ng-1)){
    for(j in (i+1):ng){
      nz.id = which(X[, i]*X[, j] != 0)
      if(length(nz.id) > 2){
        A[i, j] = cor(X[nz.id, i], X[nz.id, j])*length(nz.id)/n
        A[j, i] = A[i, j]
      }
    }
  }
  if(!is.null(colnames(X))){
    colnames(A) <- rownames(A) <- colnames(X)
  }
  return(A)
}

sLED.nonzero = function (X, Y, rho = 1000, sumabs.seq = 0.2,
          npermute = 100, useMC = FALSE, mc.cores = 1, seeds = NULL,
          verbose = TRUE, niter = 20, trace = FALSE){
  D.hat <- getDiffMatrix.nonzero(X, Y)
  pma.results <- sLEDTestStat(Dmat = D.hat, rho = rho, sumabs.seq = sumabs.seq,
                              niter = niter, trace = trace)
  Tn <- pma.results$stats
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Z <- rbind(X, Y)
  permute.results <- sLEDpermute.nonzero(Z = Z, n1 = n1, n2 = n2, rho = rho, sumabs.seq = sumabs.seq, npermute = npermute,
                                 useMC = useMC, mc.cores = mc.cores, seeds = seeds, verbose = verbose,
                                 niter = niter, trace = trace)
  pVal = rowSums(permute.results$Tn.permute > Tn)/npermute
  return(c(pma.results, permute.results, list(Tn = Tn, pVal = pVal)))
}

getDiffMatrix.nonzero = function(X, Y){
  D.hat = cor.nonzero(X) - cor.nonzero(Y)
  return(D.hat)
}

sLEDpermute.nonzero <- function(Z, n1, n2, rho=1000,sumabs.seq=0.2, npermute=100,
                                useMC=FALSE, mc.cores=1, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE) {
  ## permutation
  if (verbose) {
    cat(npermute, "permutation started:\n")
  }

  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    perm.results <- parallel::mclapply(1:npermute, sLEDOnePermute.nonzero,
                                       Z=Z, n1=n1, n2=n2, seeds=seeds,
                                       sumabs.seq=sumabs.seq,
                                       rho=rho, verbose=FALSE, niter=niter, trace=trace,
                                       mc.cores=mc.cores, mc.preschedule=TRUE)
    } else {
      ## without parallelization
      perm.results <- lapply(1:npermute, sLEDOnePermute.nonzero,
                             Z=Z, n1=n1, n2=n2, seeds=seeds,
                             sumabs.seq=sumabs.seq, rho=rho,
                             verbose=verbose, niter=niter, trace=trace)
    }

  ## extract test statistics and signs
  ntest <- length(sumabs.seq)
  Tn.permute <- matrix(NA, nrow=ntest, ncol=npermute)
  Tn.permute.sign <- matrix(NA, nrow=ntest, ncol=npermute)
  for (i in 1:npermute) {
    Tn.permute[, i] <- perm.results[[i]]$Tn.permute
    Tn.permute.sign[, i] <- perm.results[[i]]$Tn.permute.sign
  }

  if (verbose) {
    cat("permutations finished.", fill=TRUE)
  }

  return(list(Tn.permute = Tn.permute, Tn.permute.sign = Tn.permute.sign))
}

sLEDOnePermute.nonzero <- function(i, Z, n1, n2, seeds=NULL, sumabs.seq=0.2, rho=1000,
                                   verbose=TRUE, niter=20, trace=FALSE) {
  if (!is.null(seeds)){
    set.seed(seeds[i])
  }

  i.permute <- permuteIndex(n1, n2)
  D.permute <- getDiffMatrix.nonzero(X=Z[i.permute$i1, ], Y=Z[i.permute$i2, ])
  test.permute <- sLEDTestStat(Dmat=D.permute, rho=rho,
                               sumabs.seq=sumabs.seq, niter=niter, trace=trace)

  if (verbose  && (i %% 10)==0){
    cat(i, ",")
  }
  return(list(Tn.permute=test.permute$stats,
              Tn.permute.sign=test.permute$sign))
}
