#### R version 3.3.2 (2016-10-31)
#### Bioconductor 3.3 (BiocInstaller 1.22.3)
####                                         
#### by Heng Pan at Weill Cornell Medicine
####
#### This script is used to generate p values indicating patient specific epidriver inference

rm(list=ls())
library(betareg)

# load input matrix
dir <- '~/Documents/projects/CLL/git/code/' 
input.mat <- readRDS(paste0(dir, '/testdata/input-matrix.rds'))

# data transformation from [0,1] to (0,1)
scale.factor <- round(nrow(input.mat)/length(unique(input.mat$id)))
input.mat$dhcrTBeta <- (input.mat$dhcrT*(scale.factor-1)+0.5)/scale.factor

# fit beta regression model
fit.beta <- betareg(dhcrTBeta~dhcrN+pdrN+gexpN+reptime+pdrT+depthT+ncpgT, data=input.mat)

# estimate expected hypermethylation of tumor samples
input.mat$beta.response <- predict(fit.beta, data=input.mat, type='response') # mu
input.mat$beta.precision <- predict(fit.beta, data=input.mat, type='precision') # phi
shape1 <- input.mat$beta.response * input.mat$beta.precision
shape2 <- input.mat$beta.precision - shape1
ans <- data.frame(beta=input.mat$dhcrTBeta, shape1=shape1, shape2=shape2)

# evaluate if observed DHcR is significantly higher than expected DHcR
input.mat$pval <- apply(ans, 1, function(x) pbeta(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), lower.tail=F))
saveRDS(input.mat, paste0(dir, '/testdata/pval-by-gene-pt.rds'))