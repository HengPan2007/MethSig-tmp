# MethSig
Process bisulfite sequencing data and performs DNAme driver inference.

### RRBS data alignment and processing
RRBS data were aligned and processed as described in our published book chapter (Pan et al., Cancer Systems Biology, 2018).

### Calculate DHcR
Promoter (defined as ± 2kb windows centered on Refseq transcription start site) hypermethylation was measured using differentially hypermethylated cytosine ratio (DHcR), defined as the ratio of hypermethylated cytosines (HCs) to the total number of promoter CpGs profiled. HCs of each sample were defined as CpGs at which DNAme is statistically higher than the average DNAme of control samples (false discovery rate=20%, Chi-squared test). Only CpGs with read depth greater than 10 reads were included in the analysis. RRBS data of matched normal tissues were used as control samples. 

`DMC.(sample).txt` file can be used to calculate promoter DHcR, which is derived from the pipeline described in our published book chapter mentioned above. `DMC.(sample).txt` has the following columns

| Column |
| ------ |
| chr |
| pos |
| numC in control |
| numC+numT in control |
| numC in tumor |
| numC+numT in tumor |
| methylation ratio at CpG as tumor methylation / control methylation |
| Chi-square pvalue |
| adjusted pvalue |

### Calculate PDR
If all the CpGs on a specific read are methylated, or all of the CpGs on a read are unmethylated, the read is classified as concordant; otherwise it is classified as discordant. At each CpG, the PDR is equal to the number of discordant reads that cover that location divided by the total number of read that cover that location. The PDR of promoter is given by averaging the values of individual CpGs, as calculated for all CpGs within the promoter of interest with read depth greater than 10 reads and that are covered by reads that contain at least 4 CpGs.

`pdrCall_from_Bismark.py` can be used to call PDR of single CpG from bismark outputs (files starting with CpG_OB or CpG OT). The output file (`pdr.(sample).txt`) can be used to calculate promoter PDR. `pdr.(sample).txt` has the following columns

| Column |
| ------ |
| chr |
| start |
| strand |
| ConMethReadCount |
| ConUMethReadCount |
| DisReadCount |
| NAReadCount |

### Make input matrix
`makeMatrix.R` is used to make input matrix. Input files are Z-score normalzied covariates matrix (`CVMatrix-normalized.rds`) including dhcrN, pdrN, gexpN and reptime. Also, `DMC.(sample).txt` and `pdr.(sample).txt` are needed for each tumor in the cohort. Output file includes the following columns:

| Column | Description |
| ------ | ----------- |
| hugo | hugo gene symbol |
| id | tumor sample id |
| dhcrN | promoter DHcR of normal samples |
| pdrN | promoter PDR of normal samples |
| gexpN | gene expression level of normal samples |
| reptime | DNA replication time |
| pdrT | promoter PDR of tumor samples |
| depthT | promoter sequencing depth of tumor samples |
| ncpgT | number of CpGs in promoter of tumor samples |

### Patient sepcific hypermethylation inference
`pBeta.R`is used to estimate expected hypermethylation of tumor sample (expected DHcR) and evaluate if observed DHcR is significantly higher than expected DHcR.

### Tumor prevalent DNAme driver inference
`pCombine.R` is used to determine if promoter hypermethylation if overrepresented in patients (DNAme driver).

### Output table
| Column | Description |
| ------ | ----------- |
| hugo | Hugo gene symbol |
| rank | rank of each promoter based on its combined p value |
| sampleSize | number of samples with enough sequencing coverage (at least 5 CpGs with minimum 10X coverage) |
| pvalue | combined p value |
| padjust | Benjamini-Hochberg adjusted p value |

### Example
## Step 1: Prepare input matrix

`makeMatrix.R` loaded CVmatrix 'CVMatrix-normalized.rds' and built the input matrix based on an single example `SRR2069925`. Users need to prepare their own `CVMatrix-normalzied.rds`, `DMC.(sample).txt` and `pdr.(sample).txt` and replace `X` and `names.list`. Input matrix will be saved in `/testdata/input-matrix-example.rds` by default.

## Step 2: Patient specific hypermethylation inference
`pBeta.R` loaded input matrix `input-matrix.rds` and output patient specific hypermethylation inference into `/testdata/pval-by-gene-pt.rds`. Users need to prepare their own input matrix based on the instruction from Step 1.

## Step 3: Tumor prevalent DNAme driver inference
`pCombine.R` loaded patient specific hypermethylation inference matrix `pval-by-gene-pt.rds` and output tumor prevalent DNAme driver inference into `/testdata/methSig-output.txt`. Users need to prepare their own input matrix based on the instruction from Step 1.
