### RNA-seq data analysis ###

## load pacakges ##
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(pheatmap)

# read in count matrix #
Sample_count <- readRDS("Sample_count.rds")

# read in sample metadata #
sampleTable <- readRDS("sampleTable.rds")

# DEG analysis with DESeq2 #
dds <- DESeqDataSetFromMatrix(countData = Sample_count,
                              colData = sampleTable,
                              design = ~ ACC_No + Timepoint)
dds <- DESeq(dds)
res <- results(dds)
res_table <- as.data.frame(res)

# GSEA #
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
hs_msigdb_df_filter <- filter(hs_msigdb_df, hs_msigdb_df$gs_cat=="H") %>% dplyr::select("gs_name","human_gene_symbol")
msigdb_hallmark <- split(hs_msigdb_df_filter$human_gene_symbol, hs_msigdb_df_filter$gs_name)

for(i in 1:length(msigdb_hallmark)){
  msigdb_hallmark[[i]] <- unique(msigdb_hallmark[[i]])
}

res_table_filter <- filter(res_table, !is.na(res_table$stat))
res_table_filter$Gene <- rownames(res_table_filter)
Post_vs_Pre_genes <- res_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Post_vs_Pre_genes <- tibble::deframe(Post_vs_Pre_genes)

set.seed(42)
Post_vs_Pre_genes_fgseaRes <- fgsea(pathways=msigdb_hallmark, stats=Post_vs_Pre_genes, nPermSimple = 10000)

Hallmark_geneset_coverage <- data.frame(pathway = names(msigdb_hallmark), coverage = 0)
for(i in 1:length(msigdb_hallmark)){
  length = length(unique(intersect(res_table_filter$Gene, msigdb_hallmark[[i]])))
  total = length(unique(msigdb_hallmark[[i]]))
  Hallmark_geneset_coverage$coverage[which(Hallmark_geneset_coverage$pathway==names(msigdb_hallmark)[i])] = length/total
}
Hallmark_geneset_coverage$coverage <- Hallmark_geneset_coverage$coverage*100

Post_vs_Pre_genes_fgseaRes_significant <- filter(Post_vs_Pre_genes_fgseaRes, Post_vs_Pre_genes_fgseaRes$padj < 0.05)
Post_vs_Pre_genes_fgseaRes_significant <- left_join(Post_vs_Pre_genes_fgseaRes_significant, Hallmark_geneset_coverage, by="pathway")
Post_vs_Pre_genes_fgseaRes_significant$pathway <- factor(Post_vs_Pre_genes_fgseaRes_significant$pathway,
                                                         levels = Post_vs_Pre_genes_fgseaRes_significant$pathway[order(Post_vs_Pre_genes_fgseaRes_significant$NES)])


ggplot(Post_vs_Pre_genes_fgseaRes_significant, aes(x = NES, y = pathway, color = padj, size=coverage)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 6)) +  
  theme_minimal() +
  labs(x = "Pre-treatment  <--  NES  -->  Post-treatment", y = NULL,
       color = "padj", size = "Coverage (%)")