
source('scripts/impute_by_col.R')
source('scripts/normalize_by_col.R')

#pdf('derived/tables/norm_imputation_graphs.pdf', width = 7, height = 5)

# read in data
df <- read.table('extdata/210714_ms_clean.txt', sep = '\t', header= T)
ints <- as.data.frame(df[,7:30])
ints <- as.data.frame(sapply(ints, as.numeric))
ints[ints==0] <- NA

## log2, normalize and impute
ints_df <- log2(ints)
long = reshape2::melt(ints_df)
p1 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('Log2 normalized')


## subtract the median
ints_df <- normalize_by_col(ints_df) 
long = reshape2::melt(ints_df)
p2 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('With imputation + normalized')
fwrite(ints_df, 'extdata/2101714_data_clean_no_imputation.txt', sep = '\t')

## impute left shifted gaussian
ints_norm <- ints_df
ints_df <- impute_by_col(ints_df) # left-shited gaussian
long = reshape2::melt(ints_df)
p3 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('With imputation + normalized + imputed')

#print(p1)
#print(p2)
#print(p3)
graphics.off()

# load names
d_names <- fread('derived/tables/datanames_new.csv')
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')
conditions <- unique(d_names$V5)

# add columns for imputed values
for (cond in conditions){
  print(cond)
  bool_cond <- grepl(cond, colnames(ints_norm))
  ints_imputed <- ints_norm[,bool_cond]
  ints_df[[paste0('Imputed_',cond)]] <- rowSums(is.na(ints_imputed))
}

# write data
nrow(df)
nrow(ints_df)
outdf <- cbind(df[,!colnames(df) %in% d_names$newname], ints_df)
fwrite(outdf, 'extdata/2101714_data_imputed_clean.txt', sep = '\t')
