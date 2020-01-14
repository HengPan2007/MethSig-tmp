# MethSig
Process bisulfite sequencing data and performs epidriver inference.

## Getting Ready
#### Input matrix
| Covariate | Description |
| ------ | ----------- |
| dhcrN | promoter DHcR of normal samples |
| pdrN | promoter PDR of normal samples |
| gexpN | gene expression level of normal samples |
| reptime | DNA replication time |
| pdrT | promoter PDR of tumor samples |
| depthT | promoter sequencing depth of tumor samples |
| ncpgT | number of CpGs in promoter of tumor samples |

*z-score normalization is performed to all the covariates in the input matrix.

## Usage
pBeta: Estimate expected hypermethylation of tumor sample (expected DHcR) and evaluate if observed DHcR is significantly higher than expected DHcR.


pCombine: Determine if promoter hypermethylation if overrepresnted in patients (epidriver).
