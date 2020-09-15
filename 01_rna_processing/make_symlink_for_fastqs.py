# In this script, we want to create a symlink for the fastq files.
# In addition, we want to rename the fastqs:
# From: 7t_GTCCGC_L004_R1_001.fastq.gz to 7t_GTCCGC_L004_001_R1.fastq.gz
import os
import sys

original_dir = sys.argv[1]
new_dir = sys.argv[2]


for f in os.listdir(original_dir):
    if f.endswith('fastq.gz'):
        items = f.split('.')[0].split('_')
        new_filename = items[0] +  '_' + items[1] + '_' + items[2] + '_' + items[4] + '_' +  items[3] + '.fastq.gz'
        original_filepath = os.path.join(original_dir, f)
        new_filepath = os.path.join(new_dir, new_filename)
        os.symlink(original_filepath, new_filepath)
