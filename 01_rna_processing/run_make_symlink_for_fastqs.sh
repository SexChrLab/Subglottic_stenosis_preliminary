#!/bin/bash
for dir in Sample_8s1 Sample_7t Sample_5t Sample_3t Sample_1s Sample_13s Sample_9t Sample_7s Sample_6t Sample_4s2 Sample_12s Sample_6s1 Sample_4t2 Sample_15s Sample_10t Sample_8t Sample_5s1 Sample_1t Sample_3s; do

python /scratch/tphung3/SubglotticStenosis/01_rna_processing/make_symlink_for_fastqs.py /scratch/vdinu/isgs/rnaseqFastq/${dir} /scratch/tphung3/SubglotticStenosis/01_rna_processing/fastqs

done
