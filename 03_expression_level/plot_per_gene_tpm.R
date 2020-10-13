library(ggplot2)
library(ggrepel)

# Set working directory
setwd("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/03_expression_level/")

data = read.csv("TPM/all_samples_tpm.csv")
colnames(data) = c("gene", "10t", "12s", "13s", "15s", "1s", "1t", "3s", "3t", "4s2", "4t2", "5s1", "5t", "6s1", "6t", "7s", "7t", "8s1", "8t", "9t")


# genes = c("C1QBP", "CLTA", "ETFA", "LMNB2", "MRPL49", "RCN1", "SUCLG2", "U2AF1")
# genes = c("A2M", "ACTN1", "ACTR2", "APOC3", "APOH", "ARPC1B", "ARPC2", "CA1", "CAT", "CLIC1", "CNN2", "COL12A1", "CORO1C", "FGA", "FGB", "FGG", "GDI2", "HBA1", "HBB", "HBD", "IGHG1", "LTA4H", "MDH1", "NOMO1", "PGD", "PRDX2", "SERPINA3", "SLC4A1", "SPTA1", "TF", "TLN1", "TMX3", "UBA1")
# genes = c("KRT13")
# genes = c("PGD")
# genes = c("LYZ", "MPO", "AC027373.1")
# genes = c("SPRR1B", "COL1A1", "KRT16", "FAM25A", "IVL", "LSM10", "LAPTM5")
genes = c("KRT13")

for (gene in genes) {
for (index in 1:nrow(data)) {
  if (data[index,1] == gene) {
  gene_data = data[index,]
  for_plot = data.frame(tpm=c(gene_data$`12s`, gene_data$`13s`, gene_data$`15s`, gene_data$`1t`, gene_data$`3t`, gene_data$`4t2`, gene_data$`5t`, gene_data$`6t`, gene_data$`7t`, gene_data$`1s`, gene_data$`3s`, gene_data$`4s2`, gene_data$`5s1`, gene_data$`6s1`, gene_data$`7s`), labels = c(rep("noiSGS_subglottic", 3), rep("iSGS_trachea", 6), rep("iSGS_subglottic", 6)), sample_ids = c("12s", "13s", "15s", "1t", "3t", "4t2", "5t", "6t", "7t", "1s", "3s", "4s2", "5s1", "6s1", "7s"))
  
  for_plot$labels = as.character(for_plot$labels)
  for_plot$labels = factor(for_plot$labels, levels = unique(for_plot$labels))
  # pos <- position_jitter(width = 0.5, seed = 1)
  ggplot(for_plot, aes(x=labels, y=tpm, color = labels, label = sample_ids)) +
    geom_point(
      # position = pos
    ) +
    geom_label_repel(
      # position = pos
    ) +
    stat_summary(fun=mean, geom="point", shape=18,
                 size=3, color="red") +
    theme_bw() +
    labs(y="TPM", x="Status", title = gene) +
    scale_color_manual(values = c("#000000", "#56B4E9", "#E69F00")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  ggsave(paste("TPM/per_gene_tpm_plots/", as.character(data[index, 1]), ".png"), width = 7, height = 3.5, units = "in", dpi=300)
  }
}
}