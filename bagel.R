#### covariate values
# gene expression log2(FPKM+1)
fpkm <- read.table('~/Documents/Dan/info/gene-expression.txt', header=T, sep='\t')
rownames(fpkm) <- fpkm$Hugo_Symbol
fpkm$Hugo_Symbol <- NULL
tmp <- as.data.frame(t(apply(fpkm, 1, function(x) 2^x-1)))
tmp.mean <- apply(tmp, 1, mean, na.rm=T)
tmp.mean <- log2(tmp.mean+1)
tmp.mean <- (tmp.mean - mean(tmp.mean))/sd(tmp.mean)
cov.mRNA <- data.frame(Hugo=names(tmp.mean), mRNA=tmp.mean)
# promoter PDR in normal
pdr <- read.table('~/Documents/Dan/info/normal-pdr.txt', header=T, sep='\t')
rownames(pdr) <- pdr$Hugo_Symbol
pdr$Hugo_Symbol <- NULL
tmp.mean <- apply(pdr, 1, mean, na.rm=T)
tmp.mean <- (tmp.mean - mean(tmp.mean))/sd(tmp.mean)
cov.pdr <- data.frame(Hugo=names(tmp.mean), PDR=tmp.mean)
# promoter methylation
# methyl <- read.table('~/Documents/Dan/info/gene-methylation.txt', header=T, sep='\t')
# rownames(methyl) <- methyl$Hugo_Symbol
# methyl$Hugo_Symbol <- NULL
# tmp.mean <- apply(methyl, 1, mean, na.rm=T)
# tmp.mean <- (tmp.mean - mean(tmp.mean))/sd(tmp.mean)
# cov.methyl <- data.frame(Hugo=names(tmp.mean), methyl=tmp.mean)
# gene replication time
repTime <- read.table('~/Documents/Dan/info/replication-time.txt', header=T, sep='\t')
repTime <- repTime[,c('gene', 'reptime')]
colnames(repTime)[1] <- c('Hugo')
repTime <- repTime[!is.na(repTime$reptime),]
repTime$reptime <- (repTime$reptime - mean(repTime$reptime, na.rm=))/sd(repTime$reptime)
cov.reptime <- repTime
# hyper DMCs in normal samples

V <- merge(cov.mRNA, cov.pdr)
V <- merge(V, cov.reptime)
row.names(V) <- V$Hugo
V$Hugo <- NULL

#### DMC matrix
# promoter (5kb on each side)
library(GenomicRanges)
refSeq <- read.table('~/Documents/DLBCL/index/refSeq_hg19_20151111.txt', header=TRUE, sep='\t')
refSeq <- GRanges(seqnames=Rle(refSeq$chrom), strand=Rle(refSeq$strand), ranges=IRanges(start=refSeq$txStart, end=refSeq$txEnd), symbol=refSeq$name2)
seqlevels(refSeq, force=TRUE) <- paste0('chr', c(1:22, 'X', 'Y'))
proList <- promoters(refSeq, upstream=5000, downstream=5000)
proList <- unique(proList)
proList$id <- 1:length(proList)

# DMCs in normal samples
gene.DMC <- function(name, pro){
  example <- read.table(paste0(dir, '/DMC.', name, '_Normal.txt'))
  anno <- GRanges(seqnames=Rle(example$V1), ranges=IRanges(start=example$V2, end=example$V2), label=example$V11)
  index <- countOverlaps(pro, anno) >= 5
  pro <- pro[index]
  index <- as.list(findOverlaps(pro, anno))
  count.DMC <- function(id){
    ans <- table(mcols(anno[id])$'label')
    ans <- ans[c('UP', 'DOWN', '-')]
    return(ans)
  }
  cll <- as.data.frame(t(sapply(index, count.DMC)))
  cll$id <- mcols(pro)$'id'
  cll$total <- cll$UP + cll$DOWN + cll$'-'
  cll$'-' <- NULL
  return(cll)
}
matrix.clean <- function(dummy){
  index <- apply(dummy, 1, function(x) table(is.na(x))[['FALSE']])
  dummy <- dummy[index > 2,]
  dummy$id <- NULL
  tmp <- split(dummy, dummy$symbol, drop=T)
  ans <- as.data.frame(t(sapply(tmp, function(x) apply(x[,2:ncol(x)], 2, sum, na.rm=T))))
}
# dir <- '~/Documents/Dan/ERRBS/DMC/'
list <- paste('Normal', c(paste('CD19', c(2:4, 6:8, 10:12), sep='_'), paste('CD27', c(2:4, 6, 8, 10:12), sep='_'), paste('IGD', c(2, 3, 6, 8, 10:12), sep='_')), sep='_')
# name <- as.list(list)
# normal <- lapply(name, gene.DMC, pro=proList)
# saveRDS(normal, '~/Documents/Dan/ERRBS/DMC_g_normal.rds')
normal <- readRDS('~/Dropbox/DMC_g_normal.rds')
n_bkgd <- as.data.frame(mcols(proList)[,c('symbol', 'id')])
N_bkgd <- as.data.frame(mcols(proList)[,c('symbol', 'id')])
for (i in 1:length(list)){
  dmc <- normal[[i]]
  tmp <- dmc[,c('UP', 'id')]
  colnames(tmp)[1] <- list[i]
  n_bkgd <- merge(n_bkgd, tmp, all=T)
  tmp <- dmc[,c('total', 'id')]
  colnames(tmp)[1] <- list[i]
  N_bkgd <- merge(N_bkgd, tmp, all=T)
}
n_bkgd <- matrix.clean(n_bkgd)
N_bkgd <- matrix.clean(N_bkgd)
# list <- paste0('CLL_S', c('3_TP1', 11, 14, '20_TP1', '22_TP1', 23:25, 27, '32_TP1', 36, 37, '45_TP1', '48_TP1', 49, '50_TP1', 51, '54_TP1', 
#                           '56_TP1', '57_TP1', '59_TP1', 64, 65, 67, 73, 75, 77, '82_TP1', 84, 95, 97, 101, '103_TP1', 104:107, 109, 115, 116, 
#                           125, 127, 132, 134, 138:141, '149_TP1', '150_TP1', 151, 152, 154:157, 159, 161, 162, 164:166, 172, 175,179, 
#                           189, 190, 201, 202, 206, 209:212, 215:218, 220, 222, 224, 226, 227, 229, 230, 231, 234:238, 241:248, 250:253))
# name <- as.list(list)
# cll <- lapply(name, gene.DMC, pro=proList)
# cll <- readRDS('~/Documents/Dan/DMC_g_cll.rds')
# n_hyper<- as.data.frame(mcols(proList)[,c('symbol', 'id')])
# N_hyper <- as.data.frame(mcols(proList)[,c('symbol', 'id')])
# for (i in 1:length(list)){
#   dmc <- cll[[i]]
#   tmp <- dmc[,c('UP', 'id')]
#   colnames(tmp)[1] <- list[i]
#   n_hyper <- merge(n_hyper, tmp, all=T)
#   tmp <- dmc[,c('total', 'id')]
#   colnames(tmp)[1] <- list[i]
#   N_hyper <- merge(N_hyper, tmp, all=T)
# }
# n_hyper <- matrix.clean(n_hyper)
# N_hyper <- matrix.clean(N_hyper)
index <- intersect(rownames(V), rownames(n_bkgd))
# index <- intersect(index, rownames(n_hyper))
V <- V[index,]
n_bkgd <- n_bkgd[index,]
N_bkgd <- N_bkgd[index,]
# n_hyper <- n_hyper[index,]
# N_hyper <- N_hyper[index,]

#### Bagel Method
mutsig.dmc <- function(){
  
}
D <- as.matrix(dist(V, method='euclidean', diag=FALSE, upper=FALSE, p=2))
n_bkgd_g <- apply(n_bkgd, 1, sum)
N_bkgd_g <- apply(N_bkgd, 1, sum)
normal_dmc <- cbind(n_bkgd_g, N_bkgd_g)


gene <- as.list(rownames(bagel))
local.regression.bagels <- function(gene, D=D, bagel=bagel, B_max=50, Q_min=0.05){
  local.D <- D[gene,]
  local.D <- local.D[names(local.D) != gene]
  local.D <- local.D[order(local.D)]
  local.D <- local.D[1:B_max]
  local.bagel <- bagel[names(local.D),]
  gene.bagel <- as.vector(bagel[gene,])
  Q.gene.left <- function(N2, N1=gene.bagel){
    n_1 = N1[1]
    N_1 = N1[2]
    n_2 = N2[1]
    N_2 = N2[2]
    Q.left = pbinom(q=n_1, size=N_1, prob=n_2/N_2)
    Q <- 2*min(Q.left, 1-Q.left)
    return(Q)
  }
  Q <- apply(local.bagel, 1, Q.gene.left,N1=gene.bagel)
  x <- gene.bagel[1]
  X <- gene.bagel[2]
  for (i in 1:B_max){
    if (Q[i] >= Q_min) {
      x <- x + local.bagel[i,1]
      X <- X + local.bagel[i,2]
    } else {
      break
    }
  }
  return(c(x, X))
}
bagel.X <- as.data.frame(t(sapply(gene, local.regression.bagels, D=D, bagel=bagel, B_max=50, Q_min=0.05)))
x_g <- as.vector(bagel.X[,1])
names(x_g) <- rownames(bagel)
X_g <- as.vector(bagel.X[,2])
names(X_g) <- rownames(bagel)
n_bkgd_overall <- sum(n_bkgd_g)
N_bkgd_overall <- sum(N_bkgd_g)
mu_overall <- n_bkgd_overall/N_bkgd_overall
n_bkgd_p <- apply(n_bkgd, 2, sum)
N_bkgd_p <- apply(N_bkgd, 2, sum)
mu_p <- n_bkgd_p / N_bkgd_p
f_p <- mu_p / mu_overall
f_p_N <- N_bkgd_p / mean(N_bkgd_p)
mutsig.test <- function(pt, n=n_hyper, N=N_hyper, x=x_g, X=X_g, f_p=f_p, f_p_N=f_p_N){
  n.pt <- n[,pt]
  N.pt <- N[,pt]
  x.pt <- x*f_p[pt]*f_p_N[pt]
  X.pt <- X*f_p_N[pt]
  data <- cbind(n.pt, N.pt, x.pt, X.pt)
  data <- data[data[,2]>0,]
  pval <- apply(data, 1, function(x) {
    ans <- pbinom(x[1], size=x[2], prob=x[3]/x[4], lower.tail=F)
    return(2*min(ans, 1-ans))
  })
  return(pval)
}
list <- paste('Normal', c(paste('CD19', c(2:4, 6:8, 10:12), sep='_'), paste('CD27', c(2:4, 6, 8, 10:12), sep='_'), paste('IGD', c(2, 3, 6, 8, 10:12), sep='_')), sep='_')
pt_num <- mutsig.test(pt=list[1], n=n_bkgd, N=N_bkgd, x=x_g, X=X_g, f_p=f_p, f_p_N=f_p_N)
normal.ans <- data.frame(Hugo=names(pt_num), Hyper=as.numeric(pt_num < 0.05))
colnames(normal.ans)[2] <- list[1]
for (i in 2:length(list)){
  tmp <- mutsig.test(pt=list[i], n=n_bkgd, N=N_bkgd, x=x_g, X=X_g, f_p=f_p, f_p_N=f_p_N)
  tmp <- data.frame(Hugo=names(tmp), Hyper=as.numeric(tmp < 0.05))
  colnames(tmp)[2] <- list[i]
  normal.ans <- merge(normal.ans, tmp, all=T)
}
rownames(normal.ans) <- normal.ans$Hugo
normal.ans$Hugo <- NULL
normal.pt <- apply(normal.ans, 1, mean, na.rm=T)

#### bagels: cll
D <- as.matrix(dist(V, method='euclidean', diag=FALSE, upper=FALSE, p=2))
n_hyper_g <- apply(n_hyper, 1, sum)
N_hyper_g <- apply(N_hyper, 1, sum)
bagel <- cbind(n_hyper_g, N_hyper_g)
gene <- as.list(rownames(bagel))
local.regression.bagels <- function(gene, D=D, bagel=bagel, B_max=50, Q_min=0.05){
  local.D <- D[gene,]
  local.D <- local.D[names(local.D) != gene]
  local.D <- local.D[order(local.D)]
  local.D <- local.D[1:B_max]
  local.bagel <- bagel[names(local.D),]
  gene.bagel <- as.vector(bagel[gene,])
  Q.gene.left <- function(N2, N1=gene.bagel){
    n_1 = N1[1]
    N_1 = N1[2]
    n_2 = N2[1]
    N_2 = N2[2]
    Q.left = pbinom(q=n_1, size=N_1, prob=n_2/N_2)
    Q <- 2*min(Q.left, 1-Q.left)
    return(Q)
  }
  Q <- apply(local.bagel, 1, Q.gene.left,N1=gene.bagel)
  x <- gene.bagel[1]
  X <- gene.bagel[2]
  for (i in 1:B_max){
    if (Q[i] >= Q_min) {
      x <- x + local.bagel[i,1]
      X <- X + local.bagel[i,2]
    } else {
      break
    }
  }
  return(c(x, X))
}
bagel.X <- as.data.frame(t(sapply(gene, local.regression.bagels, D=D, bagel=bagel, B_max=50, Q_min=0.05)))
x_g <- as.vector(bagel.X[,1])
names(x_g) <- rownames(bagel)
X_g <- as.vector(bagel.X[,2])
names(X_g) <- rownames(bagel)
n_hyper_overall <- sum(n_hyper_g)
N_hyper_overall <- sum(N_hyper_g)
mu_overall <- n_hyper_overall/N_hyper_overall
n_hyper_p <- apply(n_hyper, 2, sum)
N_hyper_p <- apply(N_hyper, 2, sum)
mu_p <- n_hyper_p / N_hyper_p
f_p <- mu_p / mu_overall
f_p_N <- N_hyper_p / mean(N_hyper_p)
list <- paste0('CLL_S', c('3_TP1', 11, 14, '20_TP1', '22_TP1', 23:25, 27, '32_TP1', 36, 37, '45_TP1', '48_TP1', 49, '50_TP1', 51, '54_TP1', 
                          '56_TP1', '57_TP1', '59_TP1', 64, 65, 67, 73, 75, 77, '82_TP1', 84, 95, 97, 101, '103_TP1', 104:107, 109, 115, 116, 
                          125, 127, 132, 134, 138:141, '149_TP1', '150_TP1', 151, 152, 154:157, 159, 161, 162, 164:166, 172, 175,179, 
                          189, 190, 201, 202, 206, 209:212, 215:218, 220, 222, 224, 226, 227, 229, 230, 231, 234:238, 241:248, 250:253))
pval <- c()
pt_num <- mutsig.test(pt=list[1], n=n_hyper, N=N_hyper, x=x_g, X=X_g, f_p=f_p, f_p_N=f_p_N)
pval <- c(pval, pt_num)
cll.ans <- data.frame(Hugo=names(pt_num), Hyper=as.numeric(pt_num < 0.05))
colnames(cll.ans)[2] <- list[1]
for (i in 2:length(list)){
  tmp <- mutsig.test(pt=list[i], n=n_hyper, N=N_hyper, x=x_g, X=X_g, f_p=f_p, f_p_N=f_p_N)
  pval <- c(pval, tmp)
  tmp <- data.frame(Hugo=names(tmp), Hyper=as.numeric(tmp < 0.05))
  colnames(tmp)[2] <- list[i]
  cll.ans <- merge(cll.ans, tmp, all=T)
}
pval <- pval[!is.na(pval)]
qqplot(-log10(runif(length(pval))), -log10(pval), xlab='Expected', ylab='Observed', main='Bagel Type I: CLL Based')
lines(c(0,7), c(0,7), col='red')
rownames(cll.ans) <- cll.ans$Hugo
cll.ans$Hugo <- NULL
cll.pt <- data.frame(total=apply(cll.ans, 1, function(x) table(is.na(x))[['False']]))


