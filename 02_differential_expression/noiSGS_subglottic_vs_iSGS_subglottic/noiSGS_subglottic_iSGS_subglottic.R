# ---------------------------------
# noiSGS_subglottic_iSGS_subglottic
# ---------------------------------
library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)

# Set working directory
setwd('c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/02_differential_expression/')

# Reading in count-data
# Reading in expression data
# gene count information for each sample
# each column is a sample
# each row is the raw count (expression) for that gene
counts <- read.table("counts/featureCounts_noiSGS_subglottic_vs_iSGS_subglottic.tsv", sep = "\t")
colnames(counts) = c("12s", "13s", "15s", "1s", "3s", "4s2", "5s1", "6s1", "7s")
head(counts)

# Reading in gene data
# the gene file contains information about the genes
# Geneid, Chr, Start, End, Length
genes <- read.table("genes/genesID.csv", header = TRUE, sep = ",")
head(genes)

# Reading in phenotype data
# the pheno file contains information about the samples
pheno <- read.table("phenotypes/pheno_noiSGS_subglottic_iSGS_subglottic.csv", header = TRUE, sep = ",")
head(pheno)

# Create the DGEList object using the counts and genes
dge <- DGEList(counts=counts, genes=genes)
dge$samples$status <- pheno$status
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

# Normalizing gene expression distributions
dge <- calcNormFactors(dge)
dge$samples$norm.factors

cpm <- cpm(dge)
# save cpm
genes_cpm_normalized = cbind(dge$genes, cpm)
write.table(genes_cpm_normalized, "noiSGS_subglottic_vs_iSGS_subglottic/noiSGS_subglottic_vs_iSGS_subglottic_cpm_normalized.csv", row.names = F, quote = F, sep = ",")

# Unsupervised clustering of samples
status <- as.factor(c("noiSGS", "noiSGS", "noiSGS", "iSGS", "iSGS", "iSGS", "iSGS", "iSGS", "iSGS"))
lane <- as.factor(c("2", "4", "2", "3", "1", "2", "5", "2", "1"))
id <- as.factor(c("12", "13", "15", "1", "3", "4", "5", "6", "7"))

lcpm <- cpm(dge, log=TRUE)

png("c://Users/tphung3/Documents/ASU/projects/SubglotticStenosis/manuscript/figures/noiSGS_subglottic_iSGS_subglottic_mds.png", width = 10, height = 4, units = "in", res = 300)
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
contr.matrix <- makeContrasts(disease_vs_nondisease = nondisease - disease,
                              levels=colnames(design))
head(contr.matrix) # inspect the contrast matrix

v <- voom(dge, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))
```

# Obtaining a list of all genes with values for log2foldchange and pvalues
coef = 1
gene_table <- topTable(efit, coef=coef, n=Inf, adjust.method = "BH", sort.by = "logFC")
gene_table$FC = 2^gene_table$logFC #convert logFC to FC
write.table(gene_table, 'noiSGS_subglottic_vs_iSGS_subglottic/noiSGS_subglottic_vs_iSGS_subglottic_all_genes_logFC.csv.', row.names = F, quote = F, sep=',')
