# Subglottic_stenosis_preliminary

## 01_rna_processing
1. Make symlink for the fastq files
- Python script `make_symlink_for_fastqs.py`
- For command, see script `run_make_symlink_for_fastqs.sh`. 
2. Generate config for processing rna:
  ```
  python generate_json_config_rna.py
  ```
  - This script outputs `process_rna_config.json`
3. QC, Trim, Map, and featureCounts
- Snakefile `process_rna.snakefile`

## 04_measure_expression_level
- Compute TPM
  + Snakefile `measure_expression_level.snakefile`
- Concat TPM values for all samples:
  + `python concat_tpm_across_samples.py`
  + This script returns a file `all_samples_tpm.csv`
    + Each row is a gene
    + Each column is a sample in this order: 10t, 12s, 13s, 15s, 1s, 1t, 3s, 3t, 4s2, 4t2, 5s1, 5t, 6s1, 6t, 7s, 7t, 8s1, 8t, 9t
 - Subset TPM values for genes and samples 
  + `python subset_tpm.py`
  + Output 1: `ControlVsDisease_logFCgreater2_significantp_TPM.csv`
    + rows are 43 genes with log2FC greater than absolute value of 2 and unadjusted p value less than 0.05
    + columns are 6 paired control-disease samples
  + Output 2: `NonDiseaseVsDisease_logFCgreater2_significantp_TPM.csv`
    + rows are 272 genes with log2FC greater than absolute value of 2 and unadjusted p value less than 0.05
    + columns are the 3 nondisease samples and 6 disease samples

