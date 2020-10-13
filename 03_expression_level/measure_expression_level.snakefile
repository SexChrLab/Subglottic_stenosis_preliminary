import os

configfile: "/scratch/tphung3/SubglotticStenosis/01_rna_processing/process_rna_config.json"

rule all:
    input:
        expand("TPM/{sample}_XX_HISAT2_gene_featurecounts_TPM.tsv", sample=config["rna_samples_single"]),
        expand("TPM/{sample}_XX_HISAT2_gene_featurecounts_TPM.tsv", sample=config["rna_samples_multiple"])

rule measure_expression_level_tpm:
    input:
        "/scratch/tphung3/SubglotticStenosis/01_rna_processing/featureCounts/{sample}_XX_HISAT2_gene_featurecounts.tsv"
    output:
        "TPM/{sample}_XX_HISAT2_gene_featurecounts_TPM.tsv"
    params:
        script = "/home/tphung3/softwares/tanya_repos/rnaseq_analysis_scripts/measure_expression_level.py"
    shell:
        """
        python {params.script} --featureCounts_file {input} --outfile {output} --normalization tpm
        """
