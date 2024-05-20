#Packages
library(dplyr)
library(DESeq2)
library(ggplot2)

library(tidyverse)
library(apeglm)
library(IHW)
#BiocManager::install("IHW")

#Count table importing
data <- read.table("/Users/bainadandaeva/Desktop/R/DE/GSE159298_WT-VS-JmjC_RSEM_expected_count.txt", header = TRUE)
head(data)
row.names(data) = data$GeneID
data <- subset(data, select = -c(GeneID))
colnames(data) <- c("WT_1", "WT_2", "WT_3", "WT_4", "JmjC_1", "JmjC_2", "JmjC_3", "JmjC_4")
group <- c(1, 1, 1, 1, 2, 2, 2, 2)

#Create meta object
sample_names <- colnames(data)
group <- c("WT", "WT", "WT", "WT", "JmjC", "JmjC", "JmjC", "JmjC")
metadata <- data.frame(SampleName = sample_names, Group = group)
#DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Group)

#Pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


#Note on factor levels
dds$Group <- relevel(dds$Group, ref = "WT")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

summary(res)
#########
#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Group_JmjC_vs_WT", type="apeglm")

res_filtered_padj <- subset(resLFC, !is.na(padj) & padj < 0.05)

res_filtered_padj_fc <- subset(resLFC, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)



######
######
######
