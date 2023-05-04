
setwd("/projects/ps-aspenlab/demmel/HRS/saved_ouput/")
#setwd("/projects/ps-aspenlab/demmel/HRS/saved_output")
setwd("C:/Users/Admin/Documents/tui/")

# load sequence table
seqtab.nochim <- readRDS("seqtab.nochim.rds")

# CREATING ASV FASTA FILE

# store sequences and asv headers for use manipulating tables
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")

# format the headers for FASTA file
for (i in 1:dim(seqtab.nochim)[2]) {
        asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# create and write FASTA table with sequences for each ASV
asv_fasta <- c(rbind(asv_headers, asv_seqs))
saveRDS(asv_fasta, "ASV_fasta.rds")
write(asv_fasta, "ASVs.fa")

# CREATING COUNTS TABLE

# create counts table by transposing seq table and changing row names
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)

# CREATING TAXONOMY TABLE

# remove > from headers
asv_headers_tax <- sub('.', '', asv_headers)

taxa <- readRDS("taxa.rds")

# change row names from sequences to ASV numbers
rownames(taxa) <- asv_headers_tax

write.table(taxa, "ASV_taxonomy.tsv", sep = "\t", quote = F, col.names = NA)

# final output is a fasta file (and saved R object), counts table, and taxonomy table
