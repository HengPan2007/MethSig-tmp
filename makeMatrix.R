#### R version 3.3.2 (2016-10-31)
#### Bioconductor 3.3 (BiocInstaller 1.22.3)
####                                         
#### by Heng Pan at Weill Cornell Medicine
####
#### This script is used to make input matrix for pBeta.R

rm(list=ls())

# load input matrix
dir <- '~/Documents/projects/CLL/git/code/' 
X <- readRDS(paste0(dir, '/testdata/CVMatrix-normalized.rds'))

# load gene annotations
library(GenomicRanges)
gene.anno <- readRDS('~/Documents/projects/DLBCL/index/annotated.genome.hg19.rds')
tmp <- split(gene.anno, factor(mcols(gene.anno)$'exon.anno'))
gene.anno <- c(tmp$utr5, tmp$cds, tmp$utr3, tmp$intron.utr3)
geneList <- split(gene.anno, factor(mcols(gene.anno)$'symbol'))
geneList <- range(geneList)
promoterList <- promoters(geneList, upstream=2000, downstream=2000)
pro <- unlist(promoterList)

# make matrix
names.list <- as.list(c("SRR2069925"))
makeMatrix <- function(name, pro) {
  # DHcR in tumor sample
  example <- read.table(paste0(dir, '/testdata/DMC.', name, '_Normal_sum.txt'), sep='\t')
  anno <- GRanges(seqnames=Rle(example$V1), ranges=IRanges(start=example$V2, end=example$V2))
  index <- as.list(findOverlaps(pro, anno, ignore.strand=T))
  tumor <- as.data.frame(t(sapply(index, function(x) c(table(example[x, 'V11']), mean(example[x, 'V6'])))))
  colnames(tumor)[4] <- 'cov'
  tumor$total <- tumor$UP + tumor$DOWN + tumor$'-'
  tumor$hratio <- tumor$UP/tumor$total
  tumor$Hugo <- names(pro)
  tumor$'-' <- NULL
  tumor$DOWN <- NULL
  tumor <- tumor[tumor$total>=5, ]
  colnames(tumor)[1:5] <- c('hDMCtumor', 'depthT', 'ncpgT', 'dhcrT', 'hugo')
  tumor <- tumor[, c('hugo', 'depthT', 'hDMCtumor', 'ncpgT', 'dhcrT')]
  input.mat <- merge(X, tumor)
  
  # PDR in tumor samples
  data <- read.table(paste0(dir, '/testdata/pdr.',name, '.txt'), header=T, sep='\t')
  data$totalReadCount <- data$ConMethReadCount + data$ConUMethReadCount + data$DisReadCount
  data <- data[data$totalReadCount >= 10,]
  data$pdr <- data$DisReadCount / data$totalReadCount
  anno <- GRanges(seqnames=Rle(data$chr), ranges=IRanges(start=data$start, end=data$start))
  index <- as.list(findOverlaps(pro, anno, ignore.strand=T))
  pdrCal <- function(x) {
    if (length(x) < 3) {
      return(c(NA, NA))
    } else {
      return(c(length(x), mean(data[x, 'pdr'])))
    }
  }
  pdr <- as.data.frame(do.call(rbind, lapply(index, pdrCal)))
  colnames(pdr) <- c('num', 'pdrtumor')
  pdr$Hugo <- names(pro)
  pdr <- na.omit(pdr)
  pdr$num <- NULL
  colnames(pdr) <- c('pdrT', 'hugo')
  input.mat <- merge(input.mat, pdr)
  input.mat$id <- name
  return(input.mat)
}
foo <- do.call(rbind, lapply(names.list, makeMatrix, pro=pro))
foo <- foo[,c('hugo', 'id', 'dhcrN', 'pdrN', 'gexpN', 'reptime', 'pdrT', 'depthT', 'ncpgT', 'dhcrT')]
saveRDS(foo, paste0(dir, '/testdata/input-matrix-example.rds'))