

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
d_new <- data.table(
  gene_symbol = maj_gene_clean,
  #gene_symbol_last = last_gene_clean,
  gene_symbol_all = splitted_genes_pasted, 
  mappings_n = n_mappings,
  protein_ids = d$`Protein IDs`,
  peptide_counts_unique = d$`Peptide counts (unique)`
)

# Deal with LFQ
d_names <- data.frame(oldname = sort(colnames(d)[grepl('LFQ',colnames(d))]))
d_names$data_removed <- gsub('(8)?_20210630','', d_names$oldname)
d_names$id <- as.numeric(unlist(regmatches(d_names$data_removed, gregexpr("[[:digit:]]+", d_names$data_removed))))
d_names <- d_names[order(d_names$id), ]
d_names$reps <- rep(1:3, length.out = 24)
d_names$colnames <- ''

write.csv(d_names, 'derived/tables/datanames.csv', row.names = F)




sequence <- seq(1,(length(col_names)-3), by = 3)
cols <- lapply(sequence))




