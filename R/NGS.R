### NGS data analysis ###
## mutation ##
fa=list.files(path = '../../mutect_combine_new',pattern = '*.annotated.tsv')

all=c()
for( i in 1:length(fa)){
  data=read_tsv(paste0('../../mutect_combine_new/',fa[i]))
  data$tumor_name<-unique(data$tumor_name)
  if(i==1){all=rbind(all,data)}else{
    column.keep<-intersect(colnames(data),colnames(all))
    all=rbind(all[,column.keep],data[,column.keep])
  }
}

all<-as.data.frame(all)
out=all[order(all[,'chrom'],all[,'start'],all[,'end']),]

write.table(out,"allMutation_new.txt",sep="\t",quote=F,col.names=T,row.names=F)

out.filter <- out %>% filter(judgement=='KEEP') %>%
  #filter(is.na(snp138)) %>% # remove known SNP
  filter(func.knowngene%in%c('exonic',
                             'exonic;splicing',
                             'splicing') | is.na(func.knowngene)) %>% # remove mutations in non-coding regions
  filter(exonicfunc.knowngene!='synonymous SNV' | is.na(exonicfunc.knowngene)) %>% # remove synonymous SNV
  filter(esp6500siv2_all<0.005 | is.na(esp6500siv2_all)) %>%
  filter(x1kg2015aug_max<0.01 | is.na(x1kg2015aug_max)) %>%
  filter(exac_all<0.005 | is.na(exac_all)) %>%
  filter(t_alt_count>=4) %>%
  filter(tumor_f>=0.02) %>%
  filter(n_alt_count<2 | n_alt_count/(n_alt_count+n_ref_count) <0.05) %>%
  filter((n_alt_count/(n_alt_count+n_ref_count))/tumor_f < 0.2)

write.table(out.filter,"allMutation_new.filter.txt",sep="\t",quote=F,col.names=T,row.names=F)

out.filter_exo_nonsynom <- out %>% filter(judgement=='KEEP') %>%
  #filter(is.na(snp138)) %>% # remove known SNP
  filter(func.knowngene%in%c('exonic',
                             'exonic;splicing',
                             'splicing') | is.na(func.knowngene)) %>% # remove mutations in non-coding regions
  filter(exonicfunc.knowngene!='synonymous SNV' | is.na(exonicfunc.knowngene))



## indel ##
fa2=list.files(path = '../../pindel_combine_new',pattern = '*.All.annotated.tsv')

all2=c()
for( i in 1:length(fa2)){
  data=read_tsv(paste0('../../pindel_combine_new/',fa2[i]))
  if(i==1){all2=rbind(all2,data)}else{
    column.keep<-intersect(colnames(data),colnames(all2))
    all2=rbind(all2[,column.keep],data[,column.keep])
  }
}
all2<-as.data.frame(all2)
out2=all2[order(all2[,'chrom'],all2[,'start'],all2[,'end']),]
write.table(out2,"allpindel_new.txt",sep="\t",quote=F,col.names=T,row.names=F)

out.filter2 <- out2 %>%
  filter(length<50) %>%
  filter(func.knowngene%in%c('exonic',
                             'exonic;splicing',
                             'splicing') | is.na(func.knowngene)) %>% # remove mutations in non-coding regions
  filter((exonicfunc.knowngene!='nonframeshift deletion' & exonicfunc.knowngene!='nonframeshift insertion') | is.na(exonicfunc.knowngene)) %>% # remove synonymous SNV
  filter(esp6500siv2_all<0.005 | is.na(esp6500siv2_all)) %>%
  filter(x1kg2015aug_max<0.01 | is.na(x1kg2015aug_max)) %>%
  filter(exac_all<0.005 | is.na(exac_all)) %>%
  filter(t_alt_count>=4) %>%
  filter(tumor_f>=0.1) %>%
  filter(n_alt_count<2 | n_vaf <0.05) %>%
  filter(n_vaf/tumor_f < 0.2)

write.table(out.filter2,"allindel_new.filter.txt",sep="\t",quote=F,col.names=T,row.names=F)



## Oncoprint ##
mutation<-read.table('allMutation_new.filter.txt',header = T,sep='\t',quote='',stringsAsFactors = F)
mutation$exonicfunc[mutation$func!='exonic']=mutation$func[mutation$func!='exonic']
mutation.f <- mutation %>% filter(tumor_f>=0.01) %>% select(tumor_name,gene,exonicfunc)%>%unique()

indel<-read.table('allindel_new.filter.txt',header = T,sep='\t',quote='',stringsAsFactors = F)
indel<-indel[!is.na(indel$func),]
indel$exonicfunc[indel$func!='exonic']=indel$func[indel$func!='exonic']
indel.f <- indel %>% filter(length<50,tumor_f>=0.1) %>% select(tumor_name,gene,exonicfunc)%>%unique()

variant<-rbind(mutation.f,indel.f)
variant$tumor_name<-gsub('.*-.*-(.*-.*)','\\1',variant$tumor_name)

variant$exonicfunc[variant$exonicfunc=='exonic;splicing']<-'exonic'
variant$exonicfunc[variant$exonicfunc=='unknown']<-'exonic'

gene.sel<-names(table(variant$gene))[table(variant$gene)>2]

library('reshape2') 
aa=dcast(variant,gene~tumor_name,value.var = 'exonicfunc',fun.aggregate=function(x) paste(x, collapse = ";"))
rownames(aa)=aa$gene
aa=aa[,-1]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#F0F0F0', col = NA))
  },
  `frameshift deletion` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#34A02C", col = NA))
  },
  `frameshift insertion` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = 'black', col = NA))
  },
  `stopgain` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9FB5", col = NA))
  },
  `stoploss` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#73FEFF", col = NA))
  },
  `nonsynonymous SNV` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),, gp = gpar(fill = "#89419D", col = NA))
  },
  `splicing` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#C51A8A", col = NA))
  }
)

col = c('frameshift deletion' = "#34A02C", 
        'frameshift insertion'="black",
        'stopgain'="#FB9FB5",
        'stoploss'='#73FEFF',
        'nonsynonymous SNV'="#89419D",
        'splicing'="#C51A8A")


library('ComplexHeatmap')

oncoPrint(aa, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_names_gp = gpar(fontsize = 8),
          #bottom_annotation= ha,
          pct_side = "right",
          row_names_side = "left",
          show_column_names=T,
          column_title = "Landscape of Genomic Alteration",
          #column_order=surv_order,
          heatmap_legend_param = list(title = "Alterations", at = c("frameshift deletion",'frameshift insertion', "nonsynonymous SNV", 'stopgain','stoploss','splicing'),
                                      labels = c('frameshift deletion','frameshift insertion','nonsynonymous SNV','stopgain','stoploss','splicing')
          )
)