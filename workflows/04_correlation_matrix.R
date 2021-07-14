
# load data
#path <- 'extdata/2101714_data_clean_no_imputation.txt'
path <- 'extdata/2101714_data_imputed_clean.txt'
d <- fread(path, sep = '\t')

# load names
d_names <- fread('derived/tables/datanames_new.csv')
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')
mat <- d[,d_names$newname, with = F]


lst <- list()
for (c1 in colnames(mat)){
  lst[[c1]] <- list()
  for(c2 in colnames(mat)){
    cor1 <- mat[[c1]]
    cor2 <- mat[[c2]]
    #bool <- !is.na(cor1) | !is.na(cor2)
    lst[[c1]][[c2]] <- cor(cor1, cor2)
  }
}

# combine data into long format
combine <- lapply(lst, stack)
combine <- lapply(combine, function(x){colnames(x) = c('r','ind1'); return(x)})
mat_r <- do.call(rbind, lapply(names(combine), function(name){combine[[name]]$ind2 = name; return(combine[[name]])}))
mat_r$ind1 <- factor(mat_r$ind1, levels = d_names$newname)
mat_r$ind2 <- factor(mat_r$ind2, levels = d_names$newname)

ggplot(mat_r, aes(x = ind1, y = ind2, fill = r)) +
  geom_tile() + theme_bw() + xlab('') + ylab('') +
  scale_fill_gradient(low = 'white', high = 'red') +
  labs(fill = 'Pearson correlation') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Replicate correlation across experiments')
ggsave('derived/plots/correlation_matrix_all_imputed.pdf')

