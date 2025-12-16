### RNA-seq data analysis ###
### Figure 5b ###

## load pacakges ##
library(ggplot2)

### read in data ###
Bulk_on_vs_pre_fgseaRes <- read.csv("./Bulk_on_vs_pre_fgseaRes.csv")
Bulk_on_vs_pre_fgseaRes$pathway <- factor(Bulk_on_vs_pre_fgseaRes$pathway, 
                                          levels = Bulk_on_vs_pre_fgseaRes$pathway[order(Bulk_on_vs_pre_fgseaRes$NES)])

pdf("./Bulk_On_vs_Pre_GSEA.pdf",width=8,height=6)
ggplot(Post_vs_Pre_genes_fgseaRes_significant, aes(x = NES, y = pathway, color = padj, size=coverage)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +  
  scale_size(range = c(4, 6)) +  
  theme_minimal() +
  labs(x = "Pre-treatment  <-  NES  ->  On-treatment", y = NULL,
       color = "p.adj", size = "Coverage (%)")
dev.off()
