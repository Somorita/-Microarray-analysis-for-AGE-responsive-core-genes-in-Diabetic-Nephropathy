# -Microarray-analysis-for-AGE-responsive-core-genes-in-Diabetic-Nephropathy
This repository contains all code and processed data supporting the analyses reported in the manuscript

##########################################
# GSE30122 Microarray Analysis Workflow
# Author: Somorita Baishya
# Date: 2025
##########################################

# Load required libraries
library(affy)
library(GEOquery)
library(oligo)
library(Biobase)
library(tidyverse)
library(arrayQualityMetrics)
library(limma)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(msigdbr)

# -------------------------------
# 1. Load CEL files
# -------------------------------
cel_dir <- "C:/Users/shona/Documents/30122"
cel_files <- list.celfiles(cel_dir, full.names = TRUE)
cat("CEL files found:\n")
print(cel_files)

affy_batch <- ReadAffy(filenames = cel_files)

# -------------------------------
# 2. Boxplot before normalization
# -------------------------------
boxplot(affy_batch, main = "Before RMA Normalization", 
        col = "lightblue", border = "black", notch = TRUE)

# -------------------------------
# 3. RMA normalization
# -------------------------------
rma_data <- affy::rma(affy_batch)
normalized_expr <- exprs(rma_data)
normalized_expr_df <- as.data.frame(normalized_expr)

# Boxplot after normalization
boxplot(normalized_expr, main = "After RMA Normalization", 
        col = "lightblue", border = "black", notch = TRUE)

# Save normalized expression
write.csv(normalized_expr_df, 
          file = "C:/Users/shona/Documents/30122/normalized_expression.csv", 
          row.names = TRUE)

# -------------------------------
# 4. Load GEO metadata
# -------------------------------
gse <- getGEO("GSE30122", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]
metadata <- pData(gse)

# -------------------------------
# 5. Map probe IDs to gene symbols
# -------------------------------
feature_data <- fData(gse)[, c("ID", "Gene Symbol")]

expr_with_symbols <- normalized_expr_df %>%
  rownames_to_column(var = "ID") %>%
  inner_join(feature_data, by = "ID") %>%
  filter(!str_detect(`Gene Symbol`, "///")) %>%
  filter(!is.na(`Gene Symbol`) & `Gene Symbol` != "")

# Average multiple probes per gene
average_expr <- expr_with_symbols %>%
  group_by(`Gene Symbol`) %>%
  summarise(across(starts_with("GSM"), ~mean(.x, na.rm = TRUE))) %>%
  column_to_rownames("Gene Symbol")

write.csv(average_expr, 
          file = "C:/Users/shona/Documents/30122/average_expression.csv", 
          row.names = TRUE)

# -------------------------------
# 6. DEG analysis using limma
# -------------------------------
groups <- factor(metadata$characteristics_ch1.2)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(average_expr, design)
contrast_matrix <- makeContrasts(DN - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract DEGs with |log2FC| > 0 & P < 0.01
deg <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
  rownames_to_column(var = "Gene.Symbol") %>%
  filter(abs(logFC) > 0 & P.Value < 0.01)

write.csv(deg, 
          file = "C:/Users/shona/Documents/30122/DEGs_logFC0_p001.csv", 
          row.names = FALSE)

# -------------------------------
# 7. Volcano plot
# -------------------------------
ggplot(deg, aes(x = logFC, y = -log10(P.Value), 
                color = ifelse(abs(logFC) > 0 & P.Value < 0.01, "red", "black"))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot GSE30122",
       x = "Log2 Fold Change",
       y = "-log10(P-value)",
       color = "Significance") +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  theme_minimal()

# -------------------------------
# 8. GO enrichment analysis
# -------------------------------
entrez_ids <- bitr(deg$Gene.Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- entrez_ids$ENTREZID
gene_universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")

go_bp <- enrichGO(gene = gene_list, universe = gene_universe,
                  OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
dotplot(go_bp, showCategory = 15)
write.csv(go_bp, file = "C:/Users/shona/Documents/30122/GO_BP_enrichment.csv", row.names = FALSE)

# -------------------------------
# 9. KEGG enrichment
# -------------------------------
kegg_res <- enrichKEGG(gene = gene_list, organism = "hsa")
dotplot(kegg_res, showCategory = 15)
write.csv(kegg_res, file = "C:/Users/shona/Documents/30122/KEGG_enrichment.csv", row.names = FALSE)

# -------------------------------
# 10. GSEA
# -------------------------------
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
rankings <- setNames(deg$logFC, deg$Gene.Symbol)

gsea_res <- fgsea(pathways = split(hallmark_df$gene_symbol, hallmark_df$gs_name),
                  stats = rankings, minSize = 15, maxSize = 500, nperm = 10000)
head(gsea_res)
write.csv(gsea_res, file = "C:/Users/shona/Documents/30122/GSEA_results.csv", row.names = FALSE)

##########################################
# End of workflow
##########################################

