library(gplots)
library(doBy)
library(edgeR)
if (!require("doBy")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("doBy", suppressUpdates=TRUE)
  library("doBy")
}
if (!require("UpSetR")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt", suppressUpdates=TRUE)
  library(UpSetR)
}
library(reshape2)
library(RColorBrewer)
library(dplyr)

# Set working directory
setwd("../MTB1143/")

dir.create("tables")
dir.create("DE_genes_contrasts")

# Data reading and processing

#count.table <- read.table("results/featureCounts/merged_gene_counts.txt", header = T, sep = "\t", na.strings = c("","NaN"), quote=NULL, stringsAsFactors = F, dec = ".", fill=TRUE, row.names=1)
#colnames(count.table)[2:length(colnames(count.table))] <- sapply(colnames(count.table)[2:length(colnames(count.table))], function(x) unlist(strsplit(x, "_"))[1])

# Reading count table
count.table.tp1 <- read.table("mtbdata/20190125135003_mRNA_P34768Nr3.quant.tab/20190125135003_mRNA_P34768Nr3.quant.proc.tab", header = T, sep="\t")
count.table.tp2 <- read.table("mtbdata/20190125135003_mRNA_P34768Nr6.quant.tab/20190125135003_mRNA_P34768Nr6.quant.proc.tab", header = T, sep="\t")
count.table.tp1$TP1 <- count.table.tp1$unstranded
count.table.tp2$TP2 <- count.table.tp2$unstranded
count.table <- merge(x=count.table.tp1[,c("Gene","TP1")], y=count.table.tp2[,c("Gene","TP2")], by="Gene")
row.names(count.table) <- count.table$Gene
count.table$Gene <- NULL

count.table$Ensembl_ID <- row.names(count.table)

write.table(count.table, file="tables/raw_counts_gene_name.tsv", row.names = F, col.names = T, sep = "\t", quote=F)

# Make gene names table
gene_names_table <- read.table("mtbdata/20190125135003_mRNA_P34768Nr3_FPKM.tsv/20190125135003_mRNA_P34768Nr3_FPKM.tsv", header =T, sep="\t")
gene_names <- gene_names_table[, c("gene_id", "gene")]
gene_names$Ensembl_ID <- gene_names$gene_id
gene_names$gene_name <- gene_names$gene
gene_names$gene_id <- NULL
gene_names$gene <- NULL
#gene_names <- count.table[, c("Ensembl_ID", "gene_name")]

count.table$gene_name <- NULL
count.table$Ensembl_ID <- NULL



##### LOG2 FOLD CHANGE WITH edgeR
# Convert to EdgeR object DGEList 
y <- DGEList(counts=count.table, group = c(1,2))

# Filter out low expressed genes (cpm of 1 corresponds to a count of 6-7 in smallest sample)
keep <- rowSums(cpm(y)>1) >=2 
y <- y[keep,  , keep.lib.sizes=FALSE]
y$samples

# Correcting for RNA compsition (sample-specific)
y <- calcNormFactors(y)
y$samples

# Estimate dispersion (pairwise comparison between two or more groups)
# Here produces warning because there are no replicates
#y <- estimateDisp(y)

norm.counts <- cpm(y, normalized.lib.sizes = F)
norm.counts.gene <- merge(x=norm.counts, y=gene_names, by.x = "row.names", by.y = "Ensembl_ID", all.x=T)
norm.counts.gene$Ensembl_ID <- norm.counts.gene$Row.names
norm.counts.gene$Row.names <- NULL

norm.counts.gene <- norm.counts.gene %>% select("gene_name", everything())
norm.counts.gene <- norm.counts.gene %>% select("Ensembl_ID", everything())

write.table(norm.counts.gene, file = "tables/normalized_counts.tsv", row.names=F, col.names = T, sep = "\t", quote = F)

# Testing for DE genes, set dispersion 0.4 as recommended for human samples
et <- exactTest(y, dispersion=0.4, pair = c("1", "2"))


logFC_gene <- merge(x= et$table, y= gene_names, by.x = "row.names", by.y = "Ensembl_ID", all.x = T)
logFC_gene$Ensembl_ID <- logFC_gene$Row.names
logFC_gene$Row.names <- NULL

# Checks for pValue threshold and also logFC threshold
genes_dereg <- logFC_gene[which(logFC_gene$PValue <= 0.05 & (logFC_gene$logFC >= 1 | logFC_gene$logFC <= -1)),]
#genes_dereg <- logFC_gene[which(logFC_gene$logFC >= 2 | logFC_gene$logFC <= -2),]
id_dereg <- genes_dereg$Ensembl_ID

genes_dereg_table <- genes_dereg[,c("gene_name", "Ensembl_ID", "logFC", "PValue", "logCPM")]
write.table(genes_dereg_table, file = "DE_genes_contrasts/Dereg_genes_TP1_vs_TP2.tsv", row.names=F, col.names=T, sep="\t", quote = F)


