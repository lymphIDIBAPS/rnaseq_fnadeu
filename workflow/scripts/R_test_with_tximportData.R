library("tximport")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("tximportData")
library("readr")

setwd("~/RNAseq_ferran")

dir <- system.file("extdata", package = "tximportData")
list.files(dir)


samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples


files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
names(files) <- paste0("sample", 1:6)
files


tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)

sample_table_2 <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sample_table_2) <- colnames(txi.kallisto.tsv$counts)


dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sample_table_2, ~condition)

## FROM HERE TAKEN FROM pipeline_4.R
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
vsd <- assay(vsd)
# colnames(vsd) <- colData(dds)[,c("Sample")]
colnames(vsd) <- colData(dds)[,c("condition")]

# PCA
pca <- prcomp(t(vsd), scale=T)
x <- summary(pca)
pcaTable <- cbind(txi.kallisto.tsv, pca$x[,1:2])

# Plot
miniPCAtable <- pcaTable[, !colnames(pcaTable) %in% c("Seq", "Fastq1", "Fastq2", "kallisto_file")]
# meltMiniPcaTable <- melt(miniPCAtable, id.vars = c("Sample", "PC1", "PC2"))
meltMiniPcaTable <- melt(miniPCAtable, id.vars = c("sample_name", "PC1", "PC2"))

outfile = "results/pdf_test.pdf"
pdf(outfile, width = 8, height = 7, useDingbats = F)
for (var in unique(meltMiniPcaTable$variable)) { 
  p1 <- ggplot(meltMiniPcaTable[meltMiniPcaTable$variable==var,], 
               aes(x = PC1, 
                   y = PC2, 
                   label = sample_name,
                   color = value)) +
    geom_point(size=3) +
    geom_label_repel() +
    theme_classic() + 
    theme(axis.ticks = element_blank(), axis.text=element_blank()) +
    ggtitle(var)
  print(p1)
}
dev.off()