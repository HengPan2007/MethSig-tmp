#### R version 3.3.2 (2016-10-31)
#### Bioconductor 3.3 (BiocInstaller 1.22.3)
####                                         
#### by Heng Pan at Weill Cornell Medicine
####
#### This script is used to combine patient specific p values to perform epidriver inference

rm(list=ls())

# load input matrix
dir <- '~/Documents/projects/CLL/git/code/' 
foo <- readRDS(paste0(dir, '/testdata/pval-by-gene-pt.rds'))

# initiate parameters
times <- 100
K <- 10
cutoff <- 0.2

# pval combination function
pvalCombine <- function(x, K, times, cutoff=0.1) {
  # choose replicate with maximum number of covered CpGs from the same patient
  pval <- x[, c('id', 'pval', 'ncpgT')]
  tmp <- split(pval, f=pval$id, drop=T)
  pval <- sapply(tmp, function(x) x[which.max(x$ncpgT), 'pval'])
  size <- length(pval)
  p.75 <- NA
  # randomly select K (10) samples from the pool for 100 times
  if (size >= K) {
    library(metap)
    pvalSub <- function(K, pval) {
      pval.curr <- sample(pval, K)
      r <- round(cutoff*K)
      r <- ifelse(r == 0, 1, r)
      tmp <- wilkinsonp(pval.curr, r=r, alpha=0.05)
      tmp$p
    }
    pval.curr <- sapply(rep(K, times), pvalSub, pval=pval)
    p.75 <- quantile(pval.curr, probs=0.25)
  }
  return(c(size, p.75))
}
pval.by.gene <- split(foo, f=foo$hugo, drop=T)
pval <- as.data.frame(t(as.data.frame(lapply(pval.by.gene, pvalCombine, K=K, times=times, cutoff=cutoff))))
colnames(pval) <- c('sampleSize', 'pvalue')
pval$hugo <- rownames(pval)
pval <- pval[pval$sampleSize >= K,]
pval <- pval[order(pval$pvalue),]
pval$rank <- 1:nrow(pval)
pval <- pval[,c('hugo', 'rank', 'sampleSize', 'pvalue')]
pval$padjust <- p.adjust(pval$pvalue, method='BH')
write.table(pval, paste0(dir, '/testdata/methSig-output.txt'), row.names=F, sep='\t', quote=F)