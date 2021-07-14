# generate all comparisons by doing a moderated t-test
library(genoppi)

# get data
path <- 'extdata/2101714_data_imputed_clean.txt'
d <- read.csv(path, sep = '\t')
cnames = colnames(d)

# get names
d_names <- fread('derived/tables/datanames_new.csv')
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')
conditions <- unique(d_names$V5)

## contrasts
l <- list()
l$conditions <- conditions
#l <- as.list(conditions)


#l$conditions <- c(
#  'WT_DMSO', 'WT_ETO', 'WT_CPT',
#  'KO_DMSO', 'KO_ETO', 'KO_CPT',
#  'CD_DMSO', 'CD_ETO' #'CD_CPT'
#)

l$replicates <- c(1,2,3)


outlist <- list()
#pdf(paste0(path,'.pdf'), width = 14, height = 12)

for (cond1 in l$conditions) {
  outlist[[cond1]] <- list()
  for (cond2 in l$conditions) {
    
    
    if (cond1 != cond2){
      
      print(paste0(cond1,'.',cond2))
      
      # grab conditions
      bool_cond1 <- grepl(cond1, cnames) & !grepl('Imputed', cnames)
      bool_cond2 <- grepl(cond2, cnames) & !grepl('Imputed', cnames)
      
      # comparisons between conditions
      label_comparison <- paste0(cond1,'.',cond2)
      comparison <- d[,bool_cond1] - d[,bool_cond2]
      colnames(comparison) <- paste0('rep',1:3)

      # remove rows that has too many imputed values
      cond1_imputed <- d[,grepl(cond1, cnames) & grepl('Imputed', cnames)]
      cond2_imputed <- d[,grepl(cond2, cnames) & grepl('Imputed', cnames)]
      imputed <- cond1_imputed + cond2_imputed
      
      #head(comparison)
      #bool_selected <- cond1_imputed <= 1 | cond2_imputed <= 1
      bool_selected <- imputed <= 2 | cond1_imputed == 0 | cond2_imputed == 0
      #bool_selected <- cond1_imputed == 0 | cond2_imputed == 0
      #bool_selected <- T
      
      comparison <- comparison[bool_selected,]
      d_selected <- d[bool_selected, ]
      
      # calculate moderated t-test for comparison
      any(is.na(comparison))
      comparison <- calc_mod_ttest(comparison, order = F)
      fdrs <- fdrtool::fdrtool(comparison$pvalue, statistic = "pvalue", verbose = F, plot = F)
      comparison$fdrtool.qvalue <- fdrs$qval
      comparison$fdrtool.lfdr <- fdrs$lfdr
      
      # calculate coeffecient of variants for each set of replicates
      #cv_cond1 <- abs(apply(d[bool_selected,bool_cond1], 1, function(x) sd(x) / mean(x)))
      #cv_cond2 <- abs(apply(d[bool_selected,bool_cond2], 1, function(x) sd(x) / mean(x)))
      
      # get first gene-name
      #gene_df <- data.frame(gene = unlist(lapply(strsplit(d$Gene.names[bool_selected], split = ';'), function(x) sort(x)[1])))
      gene_df <- data.frame(gene = d_selected$gene_symbol, protein_ids = d_selected$protein_ids, gene_symbol_all = d_selected$gene_symbol_all)
      dnew <- cbind(gene_df, comparison)
      dnew$cond1_cell <- unlist(lapply(strsplit(cond1, split = '\\_'), function(x) x[1]))
      dnew$cond1_time <- unlist(lapply(strsplit(cond1, split = '\\_'), function(x) x[2]))
      dnew$cond2_cell <- unlist(lapply(strsplit(cond2, split = '\\_'), function(x) x[1]))
      dnew$cond2_time <- unlist(lapply(strsplit(cond2, split = '\\_'), function(x) x[2]))
      
      outlist[[cond1]][[cond2]] <- dnew
      
    }
  }
}

# combine all data
myd <- list()
for (name in names(outlist)){
  outdfs <- outlist[[name]]
  outname <- paste0(names(outdfs))
  names(outdfs) <- outname
  myd[[name]] <- outdfs
}



save(myd, file = 'extdata/myd.rda', compress = 'xz')
df <- data$fibro_240$fibro_beads %>% id_significant_proteins()
plot_volcano_basic(df) %>% make_interactive()


library(genoppi)

library(writexl)

write_xlsx(xllist,'derived/digly/210714_moderated_test.xlsx')

