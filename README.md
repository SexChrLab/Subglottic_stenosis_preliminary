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

## 02_differential_expression
### Data
- The counts files are under `02_differential_expression/counts`
- The genes id files are under `02_differential_expression/genes`
- The phenotype files are under `02_differential_expression/phenotypes`

### Codes
- iSGS_trachea_vs_iSGS_subglottic: `02_differential_expression/iSGS_trachea_vs_iSGS_subglottic/iSGS_trachea_vs_iSGS_subglottic.R`
- noiSGS_subglottic_vs_iSGS_subglottic: `02_differential_expression/noiSGS_subglottic_vs_iSGS_subglottic/noiSGS_subglottic_vs_iSGS_subglottic.R`

## 03_expression_level
- Snakefile: `03_expression_level/measure_expression_level.snakefile`
