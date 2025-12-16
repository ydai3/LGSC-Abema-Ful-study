### RPPA analysis ###
### Figure 5c ###

library(pheatmap)

### read in RPPA data ###
RPPA_z_score <- read.csv("./RPPA_z_score.csv")
RPPA_z_score <- column_to_rownames(RPPA_z_score, var="X")
RPPA_z_score <- as.matrix(RPPA_z_score)

### read in sample annotation ###
sample_anno <- read.csv("/Users/ydai6/Library/CloudStorage/OneDrive-InsideMDAnderson/Mac data/LGSOC_Abema_Ful study/code_uploading/R/RPPA/sample_anno.csv")
sample_anno <- column_to_rownames(sample_anno, var="X")
sample_anno$Timepoint <- factor(sample_anno$Timepoint, levels=c("Pre-treatment","On-treatment"))

Anno_color = list("Timepoint" = c("Pre-treatment"="#90d2c5", "On-treatment"="#fce291"))

### Plot ###
pdf("./RPPA_heatmap.pdf",width=8,height=10)
pheatmap(RPPA_z_score, fontsize_col = 12,fontsize_row = 8,
         color=colorRampPalette(c('blue','black','yellow'))(100),
         border_color = "white",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames =T,
         show_colnames =T,
         annotation_col=sample_anno,
         annotation_colors=Anno_color,
         main='Protein expression profiles')
dev.off()
