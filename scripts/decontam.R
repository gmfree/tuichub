# remember to run through the 4.0.2 version of R instead of default
library(dplyr)
packageVersion("dplyr")
library(decontam)
packageVersion("decontam")

setwd("/projects/ps-aspenlab/demmel/HRS/saved_ouput/")
#setwd("/projects/ps-aspenlab/demmel/HRS/saved_output")
setwd("C:/Users/Admin/Documents/tui/")

# load the three files generated in tables.R
counts <- read.table(file = 'ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
tax_table <- read.table(file = 'ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
asv_fasta <- readRDS("ASV_fasta.rds")

# sanity check to make sure the tables loaded right
# counts colnames should be a list of the samples, tax_table colnames should be kingdom phylum class etc...
print("Counts columns")
colnames(counts)
print("Tax columns")
colnames(tax_table)

# create a list of T/F values stating whether each sample should be treated as negatives or not
# TRUE = negatives, FALSE = samples or positives
# rep() is used to replicate the value a certain number of times - sum of all numbers in the next line should be # of total samples
#<<EDIT THE CODE BELOW BASED ON WHICH ONES ARE YOUR NEGATIVES IN ORDER>>
#vector_for_decontam <- c(rep(TRUE,5), rep(FALSE,23), rep(TRUE,4), rep(FALSE,73))
vector_for_decontam <- grepl("NEG", colnames(counts)) 


# determine contaminants and store the ASV numbers identified as contam_asvs
contam_df <- isContaminant(t(counts), neg = vector_for_decontam)
# table output tells how many of the ASVs were contaminants
print("Contaminants")
table(contam_df$contaminant)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# in out file it will display the taxonomy of each of the ASVs that were identified
# definitely check out this table to see if any contaminants look weird or wrong
tax_table[row.names(tax_table) %in% contam_asvs, ]

saveRDS(contam_asvs, "contam_asvs.rds")

# determine chloroplast, mitochondria, no phylum ASVs
chloro <- row.names(filter(tax_table, tax_table$Phylum == "Cyanobacteria" & tax_table$Order == "Chloroplast"))
mito <- row.names(filter(tax_table, tax_table$Class == "Alphaproteobacteria" & tax_table$Family == "Mitochondria"))
nobac <- row.names(filter(tax_table, tax_table$Kingdom != "Bacteria" | is.na(tax_table$Phylum)))
to.remove = c(contam_asvs, chloro, mito, nobac)
saveRDS(to.remove, "asvs_to_remove.rds")
print("Number of ASVs to remove")
length(to.remove)

# remove contaminant ASVs from the ASV fasta file
contam_indices <- which(asv_fasta %in% paste0(">", to.remove))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_decontam <- asv_fasta[- dont_want]

# remove contaminant ASVs from counts table
counts_decontam <- counts[!row.names(counts) %in% to.remove, ]

# remove contaminant ASVs from taxonomy table
tax_table_decontam <- tax_table[!row.names(tax_table) %in% to.remove, ]

# save the outputs
write(asv_fasta_decontam, "ASVs_decontam.fa")
saveRDS(asv_fasta_decontam, "ASVs_decontam.rds")
write.table(counts_decontam, "ASV_counts_decontam.tsv", sep = "\t", quote=F, col.names = NA)
write.table(tax_table_decontam, "ASV_taxonomy_decontam.tsv", sep = "\t", quote=F, col.names = NA)

# final output is a list of the contaminant ASVs (+ taxonomy), and new fasta, counts, and taxonomy files with contaminants completely removed
