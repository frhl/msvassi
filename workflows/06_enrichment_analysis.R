
library(data.table)
library(genoppi)
library(UpSetR)

# load cleaned data with comparisons
load('extdata/myd.rda')


# against beads
comp1 <- list(
  t0 = myd$fibro_0$fibro_beads, 
  t10 = myd$fibro_10$fibro_beads, 
  t20 = myd$fibro_20$fibro_beads, 
  t40 = myd$fibro_40$fibro_beads, 
  t240 = myd$fibro_240$fibro_beads 
)

# against conditions
comp2 <- list(
  t10 = myd$fibro_10$fibro_0, 
  t20 = myd$fibro_20$fibro_0, 
  t40 = myd$fibro_40$fibro_0, 
  t240 = myd$fibro_240$fibro_0 
)

# find genes that are enriched at different time points
result1 <- lapply(names(comp1), function(t){
  
  item <- comp1[[t]]
  df <- id_significant_proteins(item, fdr_cutoff = 0.1)
  df <- df[df$significant,]
  df <- df[,c('gene','logFC','pvalue','FDR')]
  df$t <- as.numeric(gsub('t','',t))
  return(df$gene)
  
})




names(result1) <- names(comp1)

#writexl::write_xlsx(result1, 'derived/tables/fibro_versus_beads.xlsx')
union_sig_genes <- unique(unlist(result1))
length(unique(unlist(result1)))





#pdf('derived/plots/upsetplot_beads.pdf',width = 10, height = 8)
upset(fromList(result1), keep.order = FALSE) 
#graphics.off()

# extract genes that are significant in at least one set
keep_genes <- unlist(lapply(names(comp2), function(t){
  
  item <- comp2[[t]]
  df <- item[item$gene %in% union_sig_genes,]
  df <- df[abs(df$logFC) > 2,]
  return(df$gene)
  
}))



fcs2 <- lapply(names(comp2), function(t){
  
  item <- comp2[[t]]
  df <- item
  df <- df[df$gene %in% keep_genes,]
  #df <- item[item$gene %in% union_sig_genes,]
  #df <- df[abs(df$logFC) > 2,]
  df <- df[,c('gene','logFC','FDR')]
  df$t <- as.numeric(gsub('t','',t))
  df <- df[order(df$gene),]
  df$cell <- 'fibro'
  return(df)
  
})


# cluster data







# long matrix
fcs_long <- as.data.frame(do.call(rbind, fcs2))
fcs_long$t <- factor(fcs_long$t, levels = unique(fcs_long$t))
fcs_long$cell <- NULL

# setup wide matrix
wide <- reshape2::dcast(fcs_long, gene ~ t, value.var = 'logFC')
wide_mat <- wide
wide_mat[is.na(wide_mat)] <- 0
rownames(wide_mat) <- wide_mat$gene
wide_mat$gene <- NULL
clust <- hclust(dist(wide_mat))


#levels <- wide_na_replaced$gene_lysine[clust$order]
#levels_gene <- unlist(lapply(strsplit(levels, split = '_'), function(x) x[1]))
#heatmap.2(wide_mat) 

fcs_long$gene <- factor(fcs_long$gene, levels = clust$labels[clust$order])
fcs_long$label <- ifelse(abs(fcs_long$logFC) > 2, '*','')

ggplot(fcs_long, aes(x=t, y=gene, fill=logFC, label = label)) +
  geom_tile() +
  geom_text() +
  scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white') + 
  theme_bw()
ggsave('derived/plots/heatmap_fdr01beads.pdf',width = 8, height = 14)



fwrite(data.frame(gene = as.character(fcs_long$gene)), file = 'derived/heatmap_fdr01beads_genes.txt')

  





















