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

#### DHcR
Promoter (defined as ± 2kb windows centered on Refseq transcription start site) hypermethylation was measured using differentially hypermethylated cytosine ratio (DHcR), defined as the ratio of hypermethylated cytosines (HCs) to the total number of promoter CpGs profiled (Fig.1a, top-left panel). HCs of each sample were defined as CpGs at which DNAme is statistically higher than the average DNAme of control samples (false discovery rate=20%, Chi-squared test). Only CpGs with read depth greater than 10 reads were included in the analysis. RRBS data of matched normal tissues were used as control samples.

## Usage
pBeta: Estimate expected hypermethylation of tumor sample (expected DHcR) and evaluate if observed DHcR is significantly higher than expected DHcR.


pCombine: Determine if promoter hypermethylation if overrepresnted in patients (epidriver).
