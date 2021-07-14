
# load libs
library(data.table)
library(ggplot2)
#library(ggrepel)

# get cleaned data
d <- fread('extdata/210714_ms_clean.txt')

# get names
d_names <- fread('derived/tables/datanames_new.csv')
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')


## all data
pdf('derived/plots/PCA_analysis1.pdf', width = 6, height = 5)

# convert to matrix
M <- as.matrix(d[,c(7:30)])
pca <- prcomp(M)
std_explained <- pca$sdev / sum(pca$sdev)

# get PC coordinates
PC <- as.data.frame(pca$rotation)
PC$newname <- rownames(PC)
PC <- merge(PC, d_names)
PC$V5 <- factor(PC$V5)

p1 <- ggplot(PC, aes(x=PC1, y=PC2, color = V5, label = reps)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle('All data')

print(p1)
#ggsave('derived/plots/')
  

#### Only fibros

# convert to matrix
M2 <- as.matrix(d[,grepl('fibro',colnames(d)), with = F])
pca2 <- prcomp(M2)
std_explained <- pca2$sdev / sum(pca2$sdev)

# get PC coordinates
PC <- as.data.frame(pca2$rotation)
PC$newname <- rownames(PC)
PC <- merge(PC, d_names)
PC$V5 <- factor(PC$V5)

p2 <- ggplot(PC, aes(x=PC1, y=PC2, color = V5, label = reps)) +
  geom_point(size = 3) +
  theme_bw() + 
  ggtitle('Fibros')
print(p2)

## Only fibros - beads

# convert to matrix
M3 <- as.matrix(d[,grepl('fibro',colnames(d)) & !grepl('beads',colnames(d)), with = F])
pca3 <- prcomp(M3)
std_explained <- pca3$sdev / sum(pca3$sdev)

# get PC coordinates
PC <- as.data.frame(pca3$rotation)
PC$newname <- rownames(PC)
PC <- merge(PC, d_names)
PC$V5 <- factor(PC$V5)

p3 <- ggplot(PC, aes(x=PC1, y=PC2, color = V5, label = reps)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle('Fibros - beads')
print(p3)


graphics.off()