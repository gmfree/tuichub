# load dada2 package

setwd("C:/Users/Admin/Documents/tui/")

library(dada2)
packageVersion('dada2')

start_time <- Sys.time() # track timing

# path to folder containing demultiplexed library sequencing files
#path <- args[1]
path = 'C:/Users/Admin/Documents/tui/trimmed/'
paste0(path) # ensure you are using the correct directory

# forward and reverse files (collect all file names- changed pattern to match)
fastqFs <- (list.files(path = path, pattern = ".*_trimmed_R1.fastq", full.names = TRUE))
fastqRs <- (list.files(path = path, pattern = ".*_trimmed_R2.fastq", full.names = TRUE))

# check same number of forward and reverse files
if(length(fastqFs) != length(fastqRs)) {stop("unequal number of forward and reverse files")}

# create lists of file paths and names for the trimmed files
trimFs <- paste0("C:/Users/Admin/Documents/tui/trimmed/", "_trimmed_R1.fastq")
trimRs <- paste0("C:/Users/Admin/Documents/tui/trimmed/", "_trimmed_R2.fastq")

# change directory to send the plots
setwd("/projects/ps-aspenlab/demmel/HRS/saved_output")
setwd("C:/Users/Admin/Documents/tui/output")

# create pdf of quality profiles for forward samples
# tried to make a loop for this but it didn't work so this is hardcoded for 105 samples
# edit lines as needed for your number of samples - four at a time looks pretty good but can definitely be condensed
pdf("qual_profiles_F.pdf")
plotQualityProfile(fastqFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("qual_profiles_R.pdf")
plotQualityProfile(fastqRs, aggregate = TRUE)
dev.off()

Sys.time() - start_time

# final output is two pdfs of quality score profiles - forward and reverse
