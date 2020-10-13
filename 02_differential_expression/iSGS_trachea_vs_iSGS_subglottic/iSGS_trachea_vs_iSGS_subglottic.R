# -------------------------------
# iSGS_trachea_vs_iSGS_subglottic
# -------------------------------
library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)

# Set working directory
setwd('c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/02_differential_expression/')

# Reading in expression data
# gene count information for each sample
# each column is a sample
# each row is the raw count (expression) for that gene
counts <- read.table("counts/featureCounts_iSGS_trachea_vs_iSGS_subglottic.tsv", sep = "\t")
colnames(counts) = c("1s", "1t", "3s", "3t", "4s2", "4t2", "5s1", "5t", "6s1", "6t", "7s", "7t")
head(counts)

# Reading in gene data
# the gene file contains information about the genes
# Geneid, Chr, Start, End, Length
genes <- read.table("genes/genesID.csv", header = TRUE, sep = ",")
head(genes)

# Reading in phenotype data
# the pheno file contains information about the samples
pheno <- read.table("phenotypes/pheno_iSGS_trachea_iSGS_subglottic.csv", header = TRUE, sep = ",")
head(pheno)

# Create the DGEList object using the counts and genes
dge <- DGEList(counts=counts, genes=genes)
dge$samples$status <- pheno$status
dge$samples$id <- pheno$id
table(dge$samples$status) # Inspecting the number (N) of samples in each group
dim(dge)

# Data pre-processing
# Transformations from the raw-scale
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE) #log=TRUE: the cpm function adds an offset to the CPM values before converting to the log2-scale.

# Removing genes that are lowly expressed
keep.exprs <- filterByExpr(dge, group=dge$samples$status)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

### C. Normalizing gene expression distributions
lcpm <- cpm(dge, log=TRUE)

dge <- calcNormFactors(dge)
dge$samples$norm.factors

cpm <- cpm(dge)
# save cpm
genes_cpm_normalized = cbind(dge$genes, cpm)
write.table(genes_cpm_normalized, "iSGS_trachea_vs_iSGS_subglottic/iSGS_trachea_vs_iSGS_subglottic_cpm_normalized.csv", row.names = F, quote = F, sep = ",")

# Unsupervised clustering of samples
status <- as.factor(c("subglottic", "trachea", "subglottic", "trachea", "subglottic", "trachea", "subglottic", "trachea", "subglottic", "trachea", "subglottic", "trachea"))
lane <- as.factor(c("3", "3", "1", "3", "2", "5", "5", "1", "2", "4", "1", "4"))
id <- as.factor(c("1", "1", "3", "3", "4", "4", "5", "5", "6", "6", "7", "7"))

lcpm <- cpm(dge, log=TRUE)

png("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/manuscript/figures/iSGS_trachea_iSGS_subglottic_mds.png", width = 10, height = 4, units = "in", res = 300)
par(mfrow=c(1,3))
col.group <- status
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set2")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
col.id <- id
levels(col.id) <-  brewer.pal(nlevels(col.id), "Set2")
col.id <- as.character(col.id)
plotMDS(lcpm, labels=status, col=col.group, dim=c(1,2), top = 100, gene.selection = "common", cex=1.4, cex.lab=1.5)
title(main="A. Sample groups Dim 1&2", cex.main=1.7)
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(1,2), top = 100, gene.selection = "common", cex=1.4, cex.lab=1.5)
title(main="B. Sequencing lanes Dim 1&2", cex.main=1.7)
plotMDS(lcpm, labels=id, col=col.id, dim=c(1,2), top = 100, gene.selection = "common", cex=1.4, cex.lab=1.5)
title(main="C. Sample IDs Dim 1&2", cex.main=1.7)
dev.off()

# Differential expression analysis
# Create a design matrix and contrasts
design <- model.matrix(~0+dge$samples$status)
colnames(design) <- gsub("dge\\$samples\\$status", "", colnames(design))
head(design)

# Contrast design for differential expression
# Defining pairwise comparisons
contr.matrix <- makeContrasts(control_vs_disease = control - disease,
                              levels=colnames(design))
head(contr.matrix) # inspect the contrast matrix

# Removing heteroscedascity from count data
v <- voom(dge, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

# Examining the number of DE genes
summary(decideTests(efit))

# Obtaining a list of all genes with values for log2foldchange and pvalues
coef = 1
gene_table <- topTable(efit, coef=coef, n=Inf, adjust.method = "BH", sort.by = "logFC")
gene_table$FC = 2^gene_table$logFC #convert logFC to FC
write.table(gene_table, 'iSGS_trachea_vs_iSGS_subglottic/iSGS_trachea_vs_iSGS_subglottic_all_genes_logFC.csv.', row.names = F, quote = F, sep=',')