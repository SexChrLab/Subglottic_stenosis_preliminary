import json
import os
from collections import defaultdict

with open('/scratch/tphung3/SubglotticStenosis/01_rna_processing/process_rna_config.json') as f:
    data = json.load(f)

rna_read_group_identifier = defaultdict(list)
rna_samples = {"rna_samples": defaultdict(list)}
rna_samples_single = {"rna_samples_single": []} #this is to store the samples with one file
rna_samples_multiple = {"rna_samples_multiple": defaultdict(list)} #this is to store the samples with multiple files

for fn in os.listdir('/scratch/tphung3/SubglotticStenosis/01_rna_processing/fastqs'):
    items = fn.split('.')[0].split('_')
    sample_id = items[0] + '_' + items[1] + '_' + items[2]
    identifier_id = items[0] + '_' + items[1] + '_' + items[2] + '_' + items[3]
    if identifier_id not in rna_samples["rna_samples"][sample_id]:
        rna_samples["rna_samples"][sample_id].append(identifier_id)
        rna_read_group_identifier["dna_read_group_identifier"].append(identifier_id)

for sample_id in rna_samples["rna_samples"]:
    if len(rna_samples["rna_samples"][sample_id]) > 1:
        for i in rna_samples["rna_samples"][sample_id]:
            rna_samples_multiple["rna_samples_multiple"][sample_id].append(i)
    else:
        for i in rna_samples["rna_samples"][sample_id]:
            rna_samples_single["rna_samples_single"].append(i)

data.update(rna_read_group_identifier)
data.update(rna_samples)
data.update(rna_samples_multiple)
data.update(rna_samples_single)

with open('/scratch/tphung3/SubglotticStenosis/01_rna_processing/process_rna_config.json', 'w') as f:
    json.dump(data, f)
