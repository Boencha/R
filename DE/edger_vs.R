data <- read.table("/Users/bainadandaeva/Desktop/R/GSE240082_gene_counts.txt", header = TRUE)
row.names(data)=data$gene_id
data <- subset(data, select = - c(gene_id))
colnames(data) <- c("H1", "H2", "N1", "N3", "N2", "H3")
group <- c(1,1,2,2,2,1)

data <- read.table("/Users/bainadandaeva/Desktop/R/DE/GSE159298_WT-VS-JmjC_RSEM_expected_count.txt", header = TRUE)
row.names(data) <- data$GeneID
data <- subset(data, select = -c(GeneID))
colnames(data) <- c("WT_1", "WT_2", "WT_3", "WT_4", "JmjC_1", "JmjC_2", "JmjC_3", "JmjC_4")
group <- c(1, 1, 1, 1, 2, 2, 2, 2)


library(edgeR)
dge <- DGEList(counts=data, genes=rownames(data), group=group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- normLibSizes(dge)


# Создание дизайна GLM
design <- model.matrix(~ group)

# Оценка дисперсии
dge <- estimateDisp(dge, design)


# Подгонка модели
fit <- glmQLFit(dge, design)

# Выполнение теста дифференциальной экспрессии
lrt <- glmQLFTest(fit)

# Получение результатов
topTags(lrt)

# Фильтрация результатов по скорректированным значениям p-уровней
sig_genes_edgeR <- lrt$table
edger_res_filtered_pval<- sig_genes_edgeR[sig_genes_edgeR$PValue < 0.05, ]
edger_res_filtered_pval_fc <- sig_genes_edgeR[sig_genes_edgeR$PValue < 0.05 & abs(sig_genes_edgeR$logFC) > 1, ]


# Построение диаграммы Венна для сигнатур
# Сначала списки генов для каждой сигнатуры, полученной с помощью DESeq2 и edgeR
genes_edgeR_padj <- rownames(edger_res_filtered_pval)
genes_DESeq2_padj <- rownames(res_filtered_padj) 

genes_edgeR_pval_fc <- rownames(edger_res_filtered_pval_fc)
genes_DESeq2_pval_fc <- rownames(res_filtered_padj_fc) 

# объекты для построения диаграммы Венна
venn_input <- list(DESeq2 = genes_DESeq2_padj, edgeR = genes_edgeR_padj)
venn_colors <- list(DESeq2 = "red", edgeR = "blue")

# Построение диаграммы Венна
library(VennDiagram)
venn.plot <- venn.diagram(
  x = venn_input,
  category.names = c("DESeq2", "edgeR"),
  filename = "venn_diagram_pval.png",
  output = TRUE,
  imagetype = "png",
  cat.col = venn_colors,
  cat.cex = 2,
  scaled = TRUE
)

venn_input_2 <- list(DESeq2 = genes_DESeq2_pval_fc, edgeR = genes_edgeR_pval_fc)
venn_colors_2 <- list(DESeq2 = "red", edgeR = "blue")

# Построение диаграммы Венна
library(VennDiagram)
venn.plot <- venn.diagram(
  x = venn_input_2,
  category.names = c("DESeq2", "edgeR"),
  filename = "venn_diagram_pval_fc.png",
  output = TRUE,
  imagetype = "png",
  cat.col = venn_colors_2,
  cat.cex = 2,
  scaled = TRUE
)
