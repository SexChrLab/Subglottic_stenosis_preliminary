import os

configfile: "process_rna_config.json"

adapter_path = "/data/storage/SAYRES/REFERENCE_GENOMES/adapters/adapter_sequence.fa"

fastq_directory = "/scratch/tphung3/SubglotticStenosis/01_rna_processing/fastqs"
fq_prefix = ["R1", "R2"]

rule all:
    input: #featureCounts
        expand("featureCounts/{sample}_XX_HISAT2_gene_featurecounts.tsv", sample=config["rna_samples_multiple"]),
        expand("featureCounts/{sample}_XX_HISAT2_gene_featurecounts.tsv", sample=config["rna_samples_single"])
    input: #merge
        expand("aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam.bai", sample=config["rna_samples_multiple"])
    input: #mapping
        expand(
			"hisat2/GRCh38.p12.genome.XXonly.{suffix}.ht2",
			suffix=[
				"1", "2", "3", "4", "5", "6", "7", "8"]),
        expand("aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.sorted.bam.bai", read_group_identifier=config["dna_read_group_identifier"]), #mapped bam index
        expand("aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.rna.sorted.stats", read_group_identifier=config["dna_read_group_identifier"]) #bam stats

    input: #after trimming qc
        "multiqc_rna_trimmed/multiqc_report.html"
    input: #Trimming
        expand("trimmed_fastqs_rna/{read_group_identifier}_trimmed_R1.fastq.gz", read_group_identifier=config["dna_read_group_identifier"]),
        expand("trimmed_fastqs_rna/{read_group_identifier}_trimmed_R2.fastq.gz", read_group_identifier=config["dna_read_group_identifier"])
    input: #before trimming qc
        "multiqc_rna/multiqc_report.html"

# ------------------
# Before trimming QC
# ------------------

rule fastqc_analysis_rna:
	input:
		os.path.join(fastq_directory, "{read_group_identifier}_{prefix}.fastq.gz")
	output:
		"fastqc_rna/{read_group_identifier}_{prefix}_fastqc.html"
	shell:
		"PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ fastqc -o fastqc_rna {input}"

rule multiqc_analysis_rna:
	input:
		expand(
			"fastqc_rna/{read_group_identifier}_{prefix}_fastqc.html",
			read_group_identifier=config["dna_read_group_identifier"], prefix=fq_prefix)
	output:
		"multiqc_rna/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f -o multiqc_rna fastqc_rna"

# --------
# Trimming
# --------
rule trim_adapters_paired_bbduk_rna:
	input:
		fq_1 = os.path.join(fastq_directory, "{read_group_identifier}_R1.fastq.gz"),
		fq_2 = os.path.join(fastq_directory, "{read_group_identifier}_R2.fastq.gz")
	output:
		out_fq_1 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R1.fastq.gz",
		out_fq_2 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R2.fastq.gz"
	params:
		adapter = adapter_path
	threads:
		2
	shell:
		"bbduk.sh -Xmx3g in1={input.fq_1} in2={input.fq_2} "
		"out1={output.out_fq_1} out2={output.out_fq_2} "
		"ref={params.adapter} ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=55 maq=20"

# ----------------
# Post trimming QC
# ----------------
rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R1.fastq.gz",
		fq2 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_rna/{read_group_identifier}_trimmed_R1_fastqc.html",
		html2 = "fastqc_trimmed_rna/{read_group_identifier}_trimmed_R2_fastqc.html"
	shell:
		"PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ fastqc -o fastqc_trimmed_rna {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed_rna:
	input:
		expand(
			"fastqc_trimmed_rna/{read_group_identifier}_trimmed_{prefix}_fastqc.html",
			read_group_identifier=config["dna_read_group_identifier"], prefix=fq_prefix)
	output:
		"multiqc_rna_trimmed/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f -o multiqc_rna_trimmed fastqc_trimmed_rna"

# Mapping with hisat2
rule hisat2_reference_index_females_dna:
	input:
		"/scratch/tphung3/PlacentaSexDiff/A_placenta/archive/batch_1_exome_processing/refs/GRCh38.p12.genome.XXonly.fa"
	output:
		expand(
			"hisat2/GRCh38.p12.genome.XXonly.{suffix}.ht2",
			suffix=[
				"1", "2", "3", "4", "5", "6", "7", "8"])
	shell:
		"hisat2-build {input} hisat2/GRCh38.p12.genome.XXonly"

rule hisat2_map_reads:
    input:
        R1 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R1.fastq.gz",
        R2 = "trimmed_fastqs_rna/{read_group_identifier}_trimmed_R2.fastq.gz"
    output:
        "aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.sorted.bam"
    params:
        threads = 8
    shell:
        "PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ hisat2 -p {params.threads} --dta "
        "-x /scratch/tphung3/SubglotticStenosis/01_rna_processing/hisat2/GRCh38.p12.genome.XXonly "
        "-1 {input.R1} -2 {input.R2} | "
        "samtools sort -O bam -o {output}"

rule index_bam_rna:
    input:
        "aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.sorted.bam"
    output:
        "aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule bam_stats_rna:
    input:
        "aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.sorted.bam"
    output:
        "aligned_bams/{read_group_identifier}.GRCh38.p12.genome.XXonly.rna.sorted.stats"
    shell:
        "samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule merge_bams_rna:
	input:
		bams = lambda wildcards: expand(
			"aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.bam",
			sample=config["rna_samples_multiple"][wildcards.sample]),
		bais = lambda wildcards: expand(
			"aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.bam.bai",
			sample=config["rna_samples_multiple"][wildcards.sample])
	output:
		"aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam"
	threads: 4
	params:
		threads = 4
	shell:
		"sambamba merge -t {params.threads} {output} {input.bams}"

rule index_merged_bams_rna:
    input:
        "aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam"
    output:
        "aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam.bai"
    shell:
        "samtools index {input}"

# featureCounts
rule xx_featurecounts_multiple:
    input:
        BAM = "aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
        GTF = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/gencode.v29.annotation.gtf"
    output:
        COUNTS = "featureCounts/{sample}_XX_HISAT2_gene_featurecounts.tsv"
    params:
        THREADS = 5
    shell:
        """
        featureCounts -T {params.THREADS} -O --primary -p -s 0 -t gene -g gene_name -a {input.GTF} -o {output.COUNTS} {input.BAM}
        """


rule xx_featurecounts_single:
    input:
        BAM = "aligned_bams/{sample}.GRCh38.p12.genome.XXonly.sorted.bam",
        GTF = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/gencode.v29.annotation.gtf"
    output:
        COUNTS = "featureCounts/{sample}_XX_HISAT2_gene_featurecounts.tsv"
    params:
        THREADS = 5
    shell:
        """
        featureCounts -T {params.THREADS} -O --primary -p -s 0 -t gene -g gene_name -a {input.GTF} -o {output.COUNTS} {input.BAM}
        """
