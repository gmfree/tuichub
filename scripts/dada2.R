# load dada2
library(dada2)
packageVersion('dada2')

setwd("/projects/ps-aspenlab/demmel/HRS/saved_ouput/")
setwd("C:/Users/Admin/Documents/tui/")

# save project directory: assumed contains directory trimmed containing trimmed reads
path = '/projects/ps-aspenlab/demmel/HRS'
path = 'C:/Users/Admin/Documents/tui/'

# creates lists of file paths and names for trimmed reads
# forward and reverse files (collect all file names- changed pattern to match)
trimFs <- (list.files(path = file.path(path, "trimmed"), pattern = ".*_trimmed_R1.fastq", full.names = TRUE))
trimRs <- (list.files(path = file.path(path, "trimmed"), pattern = ".*_trimmed_R2.fastq", full.names = TRUE))

sample.names <- gsub("_trimmed_R1.fastq", "", basename(trimFs))
#sample.names <- gsub("_trimmed_R1.fastq", "",  basename(trimFs))

# create list of file paths and names for filtered reads
filtFs <- file.path(path, "filtered", gsub("trimmed","filtered", basename(trimFs)))
filtRs <- file.path(path, "filtered", gsub("trimmed","filtered", basename(trimRs)))

# filters out reads based on quality score, truncates length, removes
# truncLen syntax is c(# bases to keep for forward reads, # bases to keep for reverse reads)
filt_out <- filterAndTrim(trimFs, filtFs, trimRs, filtRs, maxEE=c(2,2), rm.phix=TRUE, minLen=170, truncLen=c(200,175))
saveRDS(filt_out, "filt_out.rds")

# learn errors forward reads
errF <- learnErrors(filtFs, multithread = TRUE)
saveRDS(errF, "errF.rds")

# learn errors reverse reads
errR <- learnErrors(filtRs, multithread = TRUE)
saveRDS(errR, "errR.rds")

# dereplicate forward reads
derepF <- derepFastq(filtFs, verbose = TRUE)
names(derepF) <- sample.names
saveRDS(derepF, "derepF.rds")

# dereplicate reverse reads
derepR <- derepFastq(filtRs, verbose = TRUE)
names(derepR) <- sample.names
saveRDS(derepR, "derepR.rds")

# infer ASVs forward reads
dadaF <- dada(derepF, err = errF, pool = "pseudo", multithread = TRUE)
saveRDS(dadaF, "dadaF.rds")

# infer ASVs reverse reads
dadaR <- dada(derepR, err = errR, pool = "pseudo", multithread = TRUE)
saveRDS(dadaR, "dadaR.rds")

# merge amplicon pairs
merged_amplicons <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
saveRDS(merged_amplicons, "merged_amplicons.rds")

# create sequence table
seqtab <- makeSequenceTable(merged_amplicons)
saveRDS(seqtab, "seqtab.rds")

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

# tells percentage of reads retained after removing chimeras
sum(seqtab.nochim)/sum(seqtab)

# create table showing read count loss throughout the process
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names = sample.names,
dada2_input = filt_out[,1], filtered = filt_out[,2],
dada_f = sapply(dadaF, getN), dada_r = sapply(dadaR,getN),
merged = sapply(merged_amplicons, getN),
nonchim = rowSums(seqtab.nochim),
final_perc_reads_retained = round(rowSums(seqtab.nochim)/filt_out[,1]*100,1))
write.table(summary_tab, "read_count_tracking.tsv", quote = FALSE, sep = "\t", col.names = NA)

# assign taxonomy using SILVA 2019 database
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Admin/Documents/tui/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)
saveRDS(taxa, "taxa.rds")

# final output is folder of filtered reads, summary table, 10 .rds files
