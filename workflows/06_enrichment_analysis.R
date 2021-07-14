
library(data.table)
library(genoppi)

# load cleaned data with comparisons
load('extdata/myd.rda')

comp1 <- list(
  t0 = myd$fibro_0$fibro_beads, 
  t10 = myd$fibro_10$fibro_beads, 
  t20 = myd$fibro_20$fibro_beads, 
  t40 = myd$fibro_40$fibro_beads, 
  t240 = myd$fibro_240$fibro_beads 
)

# find genes that are enriched at different time points
result <- lapply(names(comp1), function(t){
  
  item <- comp1[[t]]
  df <- id_significant_proteins(item)
  df <- df[df$significant,]
  df <- df[,c('gene','logFC','pvalue','FDR')]
  df$t <- as.numeric(gsub('t','',t))
  return(df$gene)
  
})

names(result) <- names(comp1)

result

result

library(UpSetR)
pdf('derived/plots/upsetplot.pdf',width = 10, height = 8)
upset(fromList(result), keep.order = FALSE) 
graphics.off()







