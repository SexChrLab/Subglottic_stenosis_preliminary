# In this script, generate the counts.tsv file for mds. Kimberly used bash cut to extract the 7th column and paste all of the columns together.
import os
from collections import defaultdict

FILES_DIR = '/scratch/tphung3/SubglotticStenosis/04_measure_expression_level/TPM'

all_counts = []

genes_tpm = defaultdict(list)

filenames = ['10t_GTGGCC_L003', '12s_CTTGTA_L002', '13s_ATCACG_L004', '15s_TTAGGC_L002', '1s_GTAGAG_L003', '1t_CGATGT_L003', '3s_TGACCA_L001', '3t_ACAGTG_L003', '4s2_ACTTGA_L002', '4t2_GATCAG_L005', '5s1_TAGCTT_L005', '5t_GGCTAC_L001_001', '6s1_AGTCAA_L002', '6t_AGTTCC_L004_001', '7s_CAGATC_L001', '7t_GTCCGC_L004', '8s1_ATGTCA_L005_001', '8t_CCGTCC_L001', '9t_GCCAAT_L004']

for file in filenames:
    with open(os.path.join(FILES_DIR, file + '_XX_HISAT2_gene_featurecounts_TPM.tsv'), 'r') as f:
        for line in f:
            if not line.startswith('Geneid'):
                items = line.rstrip('\n').split('\t')
                genes_tpm[items[0]].append(items[8])


outfile = open('/scratch/tphung3/SubglotticStenosis/04_measure_expression_level/TPM/all_samples_tpm.csv', 'w')
header = ['Geneid', '10t', '12s', '13s', '15s', '1s', '1t', '3s', '3t', '4s2', '4t2', '5s1', '5t', '6s1', '6t', '7s', '7t', '8s1', '8t', '9t']
print (','.join(header), file=outfile)

for gene in genes_tpm:
    out = genes_tpm[gene]
    out.insert(0, gene)
    print (','.join(out), file=outfile)