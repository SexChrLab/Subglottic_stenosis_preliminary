setwd("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/02_differential_expression/")
# ----------------------------
# iSGS_trachea_iSGS_subglottic
# ----------------------------
iSGS_trachea_iSGS_subglottic = read.csv("iSGS_trachea_vs_iSGS_subglottic/iSGS_trachea_vs_iSGS_subglottic_all_genes_logFC.csv", sep=",")

# iSGS_subglottic_upregulated
iSGS_subglottic_upregulated = subset(iSGS_trachea_iSGS_subglottic, iSGS_trachea_iSGS_subglottic$FC<0.5)

write.table(iSGS_subglottic_upregulated, "iSGS_trachea_vs_iSGS_subglottic/iSGS_subglottic_upregulated.csv", quote = F, row.names = F, sep = ",")

iSGS_subglottic_upregulated_pval.filter = subset(iSGS_subglottic_upregulated, iSGS_subglottic_upregulated$P.Value<0.05)

write.table(iSGS_subglottic_upregulated_pval.filter, "iSGS_trachea_vs_iSGS_subglottic/iSGS_subglottic_upregulated_pval.filter.csv", quote = F, row.names = F, sep = ",")

# ---------------------------------
# noiSGS_subglottic_iSGS_subglottic
# ---------------------------------
noiSGS_subglottic_iSGS_subglottic = read.csv("noiSGS_subglottic_vs_iSGS_subglottic/noiSGS_subglottic_vs_iSGS_subglottic_all_genes_logFC.csv", sep=",")

# iSGS_subglottic_upregulated
iSGS_subglottic_upregulated = subset(noiSGS_subglottic_iSGS_subglottic, noiSGS_subglottic_iSGS_subglottic$FC<0.5)

write.table(iSGS_subglottic_upregulated, "noiSGS_subglottic_vs_iSGS_subglottic/iSGS_subglottic_upregulated.csv", quote = F, row.names = F, sep = ",")

iSGS_subglottic_upregulated_pval.filter = subset(iSGS_subglottic_upregulated, iSGS_subglottic_upregulated$P.Value<0.05)

write.table(iSGS_subglottic_upregulated_pval.filter, "noiSGS_subglottic_vs_iSGS_subglottic/iSGS_subglottic_upregulated_pval.filter.csv", quote = F, row.names = F, sep = ",")

# ----------------
# Find the overlap
# ----------------
iSGS_subglottic_iSGS_trachea = read.csv("iSGS_trachea_vs_iSGS_subglottic/iSGS_subglottic_upregulated_pval.filter.csv")
noiSGS_subglottic_iSGS_subglottic = read.csv("noiSGS_subglottic_vs_iSGS_subglottic/iSGS_subglottic_upregulated_pval.filter.csv")

genes_overlap = merge(iSGS_subglottic_iSGS_trachea, noiSGS_subglottic_iSGS_subglottic, by = "Geneid")