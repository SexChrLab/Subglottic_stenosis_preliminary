# genes = dict()
# with open("c://Users/tuyen/Documents/postdoc_asu/projects/SubglotticStenosis/03_differential_expression/NonDisease_vs_Disease/all_genes_logFC_NonDiseaseVsDisease.tsv", "r") as f:
#     for line in f:
#         if not line.startswith("Geneid"):
#             for line in f:
#                 items = line.rstrip("\n").split("\t")
#                 # if abs(float(items[1])) > 1:
#                 genes[items[0].split(",")[0]] = items[1]
#
# schoeff_control_up = set({"SUCLG2", "CLTA", "C1QBP", "LMNB2", "MRPL49", "ATP5B", "U2AF1", "ETFA", "RCN1"})
# schoeff_disease_up = set({"A2M", "ACTN1", "ACTR2", "APOC3", "APOH", "ARPC1B", "ARPC2", "CA1", "CAT", "CLIC1", "CNN2", "COL12A1", "CORO1C", "FGA", "FGB", "FGG", "GDI2", "HBA1", "HBB", "HBD", "IGHG1", "LTA4H", "MDH1", "NOMO1", "PGD", "PRDX2", "SERPINA3", "SLC4A1", "SPTA1", "TF", "TLN1", "TMX3", "UBA1"})
#
# # genes.intersection(schoeff_control_up)
# # genes.intersection(schoeff_disease_up)
# #
# # print(len(genes))
#
# for i in sorted(schoeff_control_up):
#     if i in genes:
#         print (i, genes[i], 2**float(genes[i]))
#     else:
#         print (i)

# ----------------------------------
# abs(log2FC) > 1 and p.value < 0.05
# ----------------------------------
genes = set()
with open("c://Users/tuyen/Documents/postdoc_asu/projects/SubglotticStenosis/03_differential_expression/NonDisease_vs_Disease/iSGS_subglottic_upregulated.tsv", "r") as f:
    for line in f:
        if not line.startswith("Geneid"):
            items = line.rstrip("\n").split("\t")
            genes.add(items[0].split(",")[0])

# schoeff_control_up = set({"SUCLG2", "CLTA", "C1QBP", "LMNB2", "MRPL49", "ATP5B", "U2AF1", "ETFA", "RCN1"})
schoeff_disease_up = set({"A2M", "ACTN1", "ACTR2", "APOC3", "APOH", "ARPC1B", "ARPC2", "CA1", "CAT", "CLIC1", "CNN2", "COL12A1", "CORO1C", "FGA", "FGB", "FGG", "GDI2", "HBA1", "HBB", "HBD", "IGHG1", "LTA4H", "MDH1", "NOMO1", "PGD", "PRDX2", "SERPINA3", "SLC4A1", "SPTA1", "TF", "TLN1", "TMX3", "UBA1"})

overlap = genes.intersection(schoeff_disease_up)
print(sorted(overlap))

