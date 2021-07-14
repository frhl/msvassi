
source('scripts/impute_by_col.R')
source('scripts/normalize_by_col.R')

pdf('derived/plots/norm_imputation_graphs.pdf', width = 7, height = 5)

# read in data
df <- read.table('extdata/210714_ms_clean.txt', sep = '\t', header= T)
ints <- as.data.frame(df[,7:30])
ints <- as.data.frame(sapply(ints, as.numeric))
ints[ints==0] <- NA

# log2, normalize and impute
ints_df <- log2(ints)
long = reshape2::melt(ints_df)
p1 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('Log2 normalized')

ints_df <- normalize_by_col(ints_df) # subtract the median
long = reshape2::melt(ints_df)
p2 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('With imputation + normalized')

ints_norm <- ints_df
ints_df <- impute_by_col(ints_df) # left-shited gaussian
long = reshape2::melt(ints_df)
p3 <- ggplot(long, aes(x = value)) +
  geom_density() +
  facet_wrap(~variable) +
  ggtitle('With imputation + normalized + imputed')

print(p1)
print(p2)
print(p3)

graphics.off()

# add columns for imputed values
for (cond in l$conditions){
  print(cond)
  bool_cond <- grepl(cond, colnames(ints))
  ints_imputed <- ints[,bool_cond]
  ints_df[[paste('Imputed',cond)]] <- rowSums(is.na(ints_imputed))
}