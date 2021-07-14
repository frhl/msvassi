# load libs
library(data.table)
library(ggplot2)

# get cleaned data
d <- fread('extdata/210714_ms_clean.txt')

# get names
d_names <- fread('derived/tables/datanames_new.csv')
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')
M <- as.matrix(d[,c(7:30)])
boxplot(M)
