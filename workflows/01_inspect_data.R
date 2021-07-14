# inspect and clean the data...
library(data.table)
source('scripts/misc.R')

# read in the data
d <- fread('extdata/proteinGroups.txt')
dim(d)

# what columns do we want?
id_cols <- c("Protein IDs", "Peptide counts (all)","Peptide counts (unique)", "Potential contaminant", "Majority protein IDs")
ids <- d[,colnames(d) %in% id_cols, with = F]

# Majority gene
maj_gene <- unlist(lapply(ids$`Majority protein IDs`, function(x) take_last(x, split = '\\|')))
maj_gene_clean <- gsub('\\_HUMAN','',maj_gene)

# get the last gene
last_gene <- unlist(lapply(ids$`Protein IDs`, function(x) take_last(x, split = '\\|')))
last_gene_clean <- gsub('\\_HUMAN','',last_gene)
ids$gene_symbol_last <- last_gene_clean

# get all genes
splitted <- strsplit(ids$`Protein IDs`,'\\;')
splitted_gene_clean <- lapply(splitted, function(x) gsub('\\_HUMAN','',unlist(lapply(x, take_last))))
n_mappings <- unlist(lapply(splitted_gene_clean, length))
splitted_genes_pasted <- unlist(lapply(splitted_gene_clean, paste, collapse = ';'))

# combine in data.frame
d_new_names <- data.table(
  gene_symbol = maj_gene_clean,
  #gene_symbol_last = last_gene_clean,
  gene_symbol_all = splitted_genes_pasted, 
  mappings_n = n_mappings,
  protein_ids = d$`Protein IDs`,
  peptide_counts_unique = d$`Peptide counts (unique)`,
  contamination = d$`Potential contaminant`
)

# Deal with LFQ names
d_names <- data.frame(oldname = sort(colnames(d)[grepl('LFQ',colnames(d))]))
d_names$data_removed <- gsub('(8)?_20210630','', d_names$oldname)
d_names$id <- as.numeric(unlist(regmatches(d_names$data_removed, gregexpr("[[:digit:]]+", d_names$data_removed))))
d_names <- d_names[order(d_names$id), ]
d_names$reps <- rep(1:3, length.out = 24)
d_names$colnames <- ''
d_names$newname <- apply(d_names[,c(3,4)], 1, paste0, collapse = '_')

# re-download them with vassi's annotation
d_names <- fread('derived/tables/datanames_new.csv')
colnames(d_names)[5] <- 'label'
d_names$newname <- apply(d_names[,c(5,4)], 1, paste0, collapse = '_')

# extract old names in order
lfq_df <- d[,d_names$oldname, with = F]
colnames(lfq_df) <- d_names$newname

# combine the data
df <- cbind(d_new_names, lfq_df)
nrow(df)

# clean up the data

## contamiantion
table(df$contamination)
df <- df[df$contamination != '+']
nrow(df)

## remove keratins
df <- df[!grepl('^KRT',df$protein_ids),] # all are removed allready

## remove ribosomal proteins
table(grepl('^RP', df$gene_symbol))
df$protein_ids[grepl('RP', df$gene_symbol)]

## at least two unique peptides
keep_peptide <- unlist(lapply(strsplit(df$peptide_counts_unique, split = ';'), function(x) any(as.numeric(x) > 1)))
table(keep_peptide)
df <- df[keep_peptide, ]


# save data
fwrite(df, 'extdata/210714_ms_clean.txt', sep = '\t')





