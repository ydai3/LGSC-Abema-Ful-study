### Figure 4 ###

library(dplyr)
library('ComplexHeatmap')

### read in mutation profile ###
mutation_profile <- read.csv("./mutation_profile.csv")
mutation_profile <- column_to_rownames(mutation_profile, var="X")

### read in sample annotation data ###
sample_anno <- read.csv("./sample_anno.csv")
sample_anno <- column_to_rownames(sample_anno, var="X")

### define sample orders ###
sample_order = colnames(mutation_profile)

### create sample annotation bar ###
anno_col <- list(
  "Timepoint" = c("Pre-treatment" = "#8dd3c7", "On-treatment" = "#fde391"),
  "Surgery_status" = c("Optimal / No residual disease" = "#a6cee3", "Suboptimal / No surgery" = "#fdbf6f"),
  "PFS_Status" = c("Long" = "#24a9e1", "Short" = "#c82127")
)

sample_top_anno <- HeatmapAnnotation(
  df = sample_anno, col = anno_col, 
  annotation_name_side = "left", show_annotation_name = TRUE
)

### define presenting format and colors for different types of alterations ###
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#F0F0F0', col = NA))
  },
  `stopgain` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9FB5", col = NA))
  },
  `nonsynonymous SNV` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),, gp = gpar(fill = "#89419D", col = NA))
  }
)

alt_col = c('stopgain'="#FB9FB5", 'nonsynonymous SNV'="#89419D")

### Plot ###
pdf("./oncoplot.pdf",width=15,height=5)
oncoPrint(mutation_profile, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = alt_col,
          row_names_gp = gpar(fontsize = 8),
          remove_empty_columns = FALSE,
          remove_empty_rows = FALSE,
          show_column_names=T,
          column_title = "Mutation landscape",
          column_order=sample_order,
          top_annotation = sample_top_anno,
          pct_side = "right",
          row_names_side = "left",
          heatmap_legend_param = list(title = "Type of mutations", at = c('nonsynonymous SNV','stopgain')
          )
)
dev.off()


