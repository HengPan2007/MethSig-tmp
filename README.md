# MethSig
Process bisulfite sequencing data and performs epidriver inference.

## Getting Ready

## Usage
pBeta: Estimate expected hypermethylation of tumor sample (expected DHcR) and evaluate if observed DHcR is significantly higher than expected DHcR.
#### Required inputs
| Covariate | Description |
| ------ | ----------- |
| `dhcrN` | run module for processing `linear` or `circ` GoT (default: `linear`) |
| `-f1/--fastqR1` | input R1.FASTQ FILE (input file can be in GZip format with .gz extension) |
| `-f2/--fastqR2` | input R2.FASTQ FILE (input file can be in GZip format with .gz extension) |
| `-c/--config` | input CONFIG FILE   (input file should be in tab-separated) |
pCombine: Determine if promoter hypermethylation if overrepresnted in patients (epidriver).
