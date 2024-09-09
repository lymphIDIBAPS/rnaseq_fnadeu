#!/usr/bin/env Rscript
# Derived from pipeline_fastq_RNAseq_3.R by Ferran Nadeu
args = commandArgs(trailingOnly=TRUE)
ensemblToOpen <- args[1]
sampleTableToOpen <- args[2]
kallistoPath <- args[3]
outfile <- args[4]

### Libraries
library("tximport")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("pepr")

### Metadata and input files

# ensembl transcripts to gene table
ensembl <- read.table(ensemblToOpen, header = T, sep = "\t", stringsAsFactors = F)
tx2gene <- ensembl[,c(4,2)]

# metadata
project_PEP = Project(file = sampleTableToOpen)
project_PEP_table = sampleTable(project_PEP)
#sampleTable <- read.table(sampleTableToOpen, sep = "\t", header = T, stringsAsFactors = F)

# kallisto files
# I DONT SEE THE USE FOR THIS CODE
kallistoFilesToOpen <- list.files(kallistoPath, "_abundance.tsv", full.names = T)
names(kallistoFilesToOpen) <- sapply(kallistoFilesToOpen, function(x) gsub("_abundance.tsv", "", strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]))
sampleTable$kallisto_file <- kallistoFilesToOpen[match(sampleTable$Sample, names(kallistoFilesToOpen))]

### Import, normalized counts and calculate PCA

# Import using tximport package
txi <- tximport(files = sampleTable$kallisto_file, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)

# DESeq
dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~1)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
vsd <- assay(vsd)
colnames(vsd) <- colData(dds)[,c("Sample")]

# PCA
pca <- prcomp(t(vsd), scale=T)
x <- summary(pca)
pcaTable <- cbind(sampleTable, pca$x[,1:2])

# Plot
miniPCAtable <- pcaTable[, !colnames(pcaTable) %in% c("Seq", "Fastq1", "Fastq2", "kallisto_file")]
meltMiniPcaTable <- melt(miniPCAtable, id.vars = c("Sample", "PC1", "PC2"))

pdf(outfile, width = 8, height = 7, useDingbats = F)
for (var in unique(meltMiniPcaTable$variable)) { 
  p1 <- ggplot(meltMiniPcaTable[meltMiniPcaTable$variable==var,], 
               aes(x = PC1, 
                   y = PC2, 
                   label = Sample,
                   color = value)) +
    geom_point(size=3) +
    geom_label_repel() +
    theme_classic() + 
    theme(axis.ticks = element_blank(), axis.text=element_blank()) +
    ggtitle(var)
  print(p1)
}
dev.off()