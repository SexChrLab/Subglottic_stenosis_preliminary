# In this script, generate the counts.tsv file for mds. Kimberly used bash cut to extract the 7th column and paste all of the columns together.
import os

FILES_DIR = '/scratch/tphung3/SubglotticStenosis/01_rna_processing/featureCounts/'

all_counts = []

# filenames = ['10t_GTGGCC_L003', '12s_CTTGTA_L002', '13s_ATCACG_L004', '15s_TTAGGC_L002', '1s_GTAGAG_L003', '1t_CGATGT_L003', '3s_TGACCA_L001', '3t_ACAGTG_L003', '4s2_ACTTGA_L002', '4t2_GATCAG_L005', '5s1_TAGCTT_L005', '5t_GGCTAC_L001_001', '6s1_AGTCAA_L002', '6t_AGTTCC_L004_001', '7s_CAGATC_L001', '7t_GTCCGC_L004', '8s1_ATGTCA_L005_001', '8t_CCGTCC_L001', '9t_GCCAAT_L004']

filenames = ['12s_CTTGTA_L002', '13s_ATCACG_L004', '15s_TTAGGC_L002', '1s_GTAGAG_L003', '3s_TGACCA_L001', '4s2_ACTTGA_L002', '5s1_TAGCTT_L005', '6s1_AGTCAA_L002', '7s_CAGATC_L001']

for file in filenames:
    with open(os.path.join(FILES_DIR, file + '_XX_HISAT2_gene_featurecounts.tsv'), 'r') as f:
        count = []
        for line in f:
            if not line.startswith('#'):
                if not line.startswith('Geneid'):
                    count.append(line.rstrip('\n').split('\t')[6])
        all_counts.append(count)

n_genes = len(all_counts[0])
print(n_genes)

outfile = open('/scratch/tphung3/SubglotticStenosis/02_mds/featureCounts_nondisease_vs_disease.tsv', 'w')

for i in range(n_genes):
    out = []
    for count in all_counts:
        out.append(count[i])
    print('\t'.join(out), file=outfile)


