# In this script, we want to subset the file all_samples_tpm.tsv to contain:
# 1. Genes where log2FC > abs(2) for 6 paired control-disease samples (genes in the file ControlVsDisease_logFCgreater2_significantp.txt)
# 2. Genes where log2FC > abs(2) for 3 nondisease and 6 disease (genes in the file NonDiseaseVsDisease_logFCgreater2_significantp.txt)

import pandas as pd
data = pd.read_csv('/scratch/tphung3/SubglotticStenosis/04_measure_expression_level/TPM/all_samples_tpm.csv')

# 6 paired control-disease sample
control_disease_gene_list = []
with open('/scratch/tphung3/SubglotticStenosis/03_differential_expression/ControlVsDisease_logFCgreater2_significantp.txt', 'r') as f:
    for line in f:
        control_disease_gene_list.append(line.rstrip('\n'))

data_control_disease = data[data['Geneid'].isin(control_disease_gene_list)]

control_disease_samples = ['Geneid', '1s', '1t', '3s', '3t', '4s2', '4t2', '5s1', '5t', '6s1', '6t', '7s', '7t']
data_control_disease_samples = data_control_disease[control_disease_samples]
print (data_control_disease_samples.shape[0])
print (data_control_disease_samples.shape[1])
data_control_disease_samples.to_csv('/scratch/tphung3/SubglotticStenosis/04_measure_expression_level/TPM/ControlVsDisease_logFCgreater2_significantp_TPM.csv', index=False)

# 3 NonDisease + 6 Disease
nondisease_disease_gene_list = []
with open('/scratch/tphung3/SubglotticStenosis/03_differential_expression/NonDiseaseVsDisease_logFCgreater2_significantp.txt', 'r') as f:
    for line in f:
        nondisease_disease_gene_list.append(line.rstrip('\n'))

data_nondisease_disease = data[data['Geneid'].isin(nondisease_disease_gene_list)]

nondisease_disease_samples = ['Geneid', '12s', '13s', '15s', '1s', '3s', '4s2', '5s1', '6s1', '7s']
data_nondisease_disease_samples = data_nondisease_disease[nondisease_disease_samples]
print (data_nondisease_disease_samples.shape[0])
print (data_nondisease_disease_samples.shape[1])
data_nondisease_disease_samples.to_csv('/scratch/tphung3/SubglotticStenosis/04_measure_expression_level/TPM/NonDiseaseVsDisease_logFCgreater2_significantp_TPM.csv', index=False)

