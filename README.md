# MethSig
Process bisulfite sequencing data and performs epidriver inference.

## Getting Ready

#### RRBS data alignment and processing
RRBS data were aligned and processed as described in our published book chapter (Pan et al., Cancer Systems Biology, 2018).

#### DHcR
Promoter (defined as ± 2kb windows centered on Refseq transcription start site) hypermethylation was measured using differentially hypermethylated cytosine ratio (DHcR), defined as the ratio of hypermethylated cytosines (HCs) to the total number of promoter CpGs profiled. HCs of each sample were defined as CpGs at which DNAme is statistically higher than the average DNAme of control samples (false discovery rate=20%, Chi-squared test). Only CpGs with read depth greater than 10 reads were included in the analysis. RRBS data of matched normal tissues were used as control samples. 

DMC.(sample).txt file can be used to calculate DHcR, which is derived from the pipeline described in our published book chapter mentioned above. DMC.(sample).txt has the following columns

| Column |
| ------ |
| chr |
| pos |
| numC in control |
| numC+numT in control |
| numC in tumor |
| numC+numT in tumor |
| methylation ratio at CpG as methylation tumor / methylation control |
| Chi-square pvalue |
| adjusted pvalue |

#### PDR
If all the CpGs on a specific read are methylated, or all of the CpGs on a read are unmethylated, the read is classified as concordant; otherwise it is classified as discordant. At each CpG, the PDR is equal to the number of discordant reads that cover that location divided by the total number of read that cover that location. The PDR of promoter is given by averaging the values of individual CpGs, as calculated for all CpGs within the promoter of interest with read depth greater than 10 reads and that are covered by reads that contain at least 4 CpGs.

#### Make input matrix
makeMatrix.R is used to make input matrix.
Input: Z-score normalzied covariates matrix including 

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


## Usage
pBeta: Estimate expected hypermethylation of tumor sample (expected DHcR) and evaluate if observed DHcR is significantly higher than expected DHcR.

pCombine: Determine if promoter hypermethylation if overrepresnted in patients (epidriver).

## Output table

| Column | Description |
| ------ | ----------- |
| hugo | Hugo gene symbol |
| rank | rank of each promoter based on its combined p value |
| sampleSize | number of samples with enough sequencing coverage (at least 5 CpGs with minimum 10X coverage |
| pvalue | combined p value |
| padjust | Benjamini-Hochberg adjusted p value |
