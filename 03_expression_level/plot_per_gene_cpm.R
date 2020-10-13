library(ggplot2)
data = read.table('c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/02_mds/pheno_six_disease_control_pairs.csv', sep = ',', header = T)

setwd('c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/02_mds/')

counts <- read.table("featureCounts_six_samples.tsv", sep = "\t")
colnames(counts) = c("1s", "1t", "3s", "3t", "4s2", "4t2", "5s1", "5t", "6s1", "6t", "7s", "7t")
head(counts)

genes <- read.table("genesID.csv", header = TRUE, sep = ",")
head(genes)

pheno <- read.table("pheno_six_disease_control_pairs.csv", header = TRUE, sep = ",")
head(pheno)

dge <- DGEList(counts=counts, genes=genes)
dge$samples$status <- pheno$status
dge$samples$id <- pheno$id
table(dge$samples$status) # Inspecting the number (N) of samples in each group
dim(dge)

cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE) #log=TRUE: the cpm function adds an offset to the CPM values before converting to the log2-scale.

keep.exprs <- filterByExpr(dge, group=dge$samples$status)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

# lcpm <- cpm(dge, log=TRUE)
cpm <- cpm(dge)
gene_cpm= cbind(dge$genes, cpm)

# genes = c("ARPC1B", "PRDX2", "FGG", "UBA1", "PGD", "HBB", "COL12A1")
genes = c("KRT13")

for (gene in genes) {
  for (index in 1:nrow(gene_cpm)) {
    if (gene_cpm[index,1] == gene) {
      gene_data = gene_cpm[index,]
      for_plot = data.frame(tpm=c(gene_data$`1s`, gene_data$`3s`, gene_data$`4s2`, gene_data$`5s1`, gene_data$`6s1`, gene_data$`7s`, gene_data$`1t`, gene_data$`3t`, gene_data$`4t2`, gene_data$`5t`, gene_data$`6t`, gene_data$`7t`), labels = c(rep("iSGS_subglottic", 6), rep("iSGS_trachea", 6)), sample_ids = c("1s", "3s", "4s2", "5s1", "6s1", "7s", "1t", "3t", "4t2", "5t", "6t", "7t"))
      
      for_plot$labels = as.character(for_plot$labels)
      for_plot$labels = factor(for_plot$labels, levels = unique(for_plot$labels))
      pos <- position_jitter(width = 0.5, seed = 1)
      ggplot(for_plot, aes(x=labels, y=tpm, color = labels, label = sample_ids)) +
        geom_point() +
        geom_label_repel() +
        stat_summary(fun=mean, geom="point", shape=18,
                     size=3, color="red") +
        theme_bw() +
        labs(y="CPM", x="Status") +
        scale_color_manual(values = c("#000000", "#E69F00")) +
        theme(legend.position = "none")
      ggsave(paste("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/04_measure_expression_level/CPM/iSGS_subglottic_iSGS_trachea/", as.character(gene_cpm[index, 1]), ".png"), width = 7, height = 3.5, units = "in")
    }
  }
}

#
data = read.table('c://Users/tuyen/Documents/postdoc_asu/projects/SubglotticStenosis/02_mds/pheno_nondisease_vs_disease.csv', sep = ',', header = T)
setwd('c://Users/tuyen/Documents/postdoc_asu/projects/SubglotticStenosis/02_mds/')

counts <- read.table("featureCounts_nondisease_vs_disease.tsv", sep = "\t")
colnames(counts) = c("12s", "13s", "15s", "1s", "3s", "4s2", "5s1", "6s1", "7s")
head(counts)

genes <- read.table("genesID.csv", header = TRUE, sep = ",")
head(genes)

pheno <- read.table("pheno_nondisease_vs_disease.csv", header = TRUE, sep = ",")
head(pheno)

dge <- DGEList(counts=counts, genes=genes)
dge$samples$status <- pheno$status
table(dge$samples$status) # Inspecting the number (N) of samples in each group
dim(dge)

cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE) #log=TRUE: the cpm function adds an offset to the CPM values before converting to the log2-scale.

keep.exprs <- filterByExpr(dge, group=dge$samples$status)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

cpm <- cpm(dge)

gene_cpm= cbind(dge$genes, cpm)

# genes = c('APOH', 'ARPC1B', 'FGB', 'HBB', 'IGHG1', 'SLC4A1', 'UBA1')
genes = c("KRT13")
for (gene in genes) {
  for (index in 1:nrow(gene_cpm)) {
    if (gene_cpm[index,1] == gene) {
      gene_data = gene_cpm[index,]
      for_plot = data.frame(tpm=c(gene_data$`1s`, gene_data$`3s`, gene_data$`4s2`, gene_data$`5s1`, gene_data$`6s1`, gene_data$`7s`, gene_data$`12s`, gene_data$`13s`, gene_data$`15s`), labels = c(rep("iSGS_subglottic", 6), rep("noiSGS_subglottic", 3)), sample_ids = c("1s", "3s", "4s2", "5s1", "6s1", "7s", "12s", "13s", "15s"))

      for_plot$labels = as.character(for_plot$labels)
      for_plot$labels = factor(for_plot$labels, levels = unique(for_plot$labels))
      pos <- position_jitter(width = 0.5, seed = 1)
      ggplot(for_plot, aes(x=labels, y=tpm, color = labels, label = sample_ids)) +
        geom_point() +
        geom_label_repel() +
        stat_summary(fun.y=mean, geom="point", shape=18,
                     size=3, color="red") +
        theme_bw() +
        labs(y="CPM", x="Status") +
        scale_color_manual(values = c("#000000", "#56B4E9")) +
        theme(legend.position = "none")
      ggsave(paste("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/04_measure_expression_level/CPM/noiSGS_subglottic_iSGS_subglottic/", as.character(gene_cpm[index, 1]), ".png"), width = 7, height = 3.5, units = "in")
    }
  }
}