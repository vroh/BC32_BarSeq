##################################################
#                                                #
#          Barcode Processing v.0.3              #
#                                                #
#               Roh et al. 2017                  #
#                                                #
#   Author: roh(dot)vincent(at)gmail(dot)com     #
#                                                #
##################################################

# This script processes a BC32 barcode sequencing run (Using R1 reads only).
# Don't forget to set your working directory to point at the directory containing sequencing data and the sampname table!

# Depends on the following packages (make sure they are installed)
library(ShortRead)
library(Biostrings)
library(dplyr)
library(stringdist)

# For indels, an updated version of Biostrings' vmatchPattern is required (described in functions.R)
source('functions.R')

# Define variables
pat <- 'CTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNN' # matching pattern (BC32)
bb_mis <- 1 # number of mismatch allowed in barcode backbone sequence
indels <- 0 # total edit distance resulting from indels that are tolerated for barcode matching (barcodes with indels increase computation time)
threshold <- 3 # threshold for Hamming distance in step pooling similar sequences

# Collect info from user
sampname <- read.delim('sampname.txt', # sampname provides a list of sample names, matching sequencing files and multiplexing index
                       header = F,
                       stringsAsFactor = F)
colnames(sampname) <- c('sample', 'file', 'index')
seqfiles <- sampname$file

# 1. Loop through samples and generate summaries
################################################

for (i in 1:length(seqfiles)) {
  
  # 1.a. Generate summaries for each sample
  #########################################
  
  # Read fastq file with ShortRead package
  sequences <- readFastq(seqfiles[i])@sread
  name <- sampname$sample[i]
  message(paste('Processing sample: ', name))
  
  # Filter for sequences with specific multiplexing index (only keep sequences with the correct index, beware of 6 and 8-mer indexes!)
  # This steps avoids misassignment of samples at de-multiplexing 
  idx <- paste0('GTCAC', sampname$index[i], 'ATCTC')
  # Index matching (strict!)
  idx_match <- vmatchPattern(idx, sequences, 0)
  # Subset sequences (and compute filtering ratio)
  init_len <- length(sequences)
  sequences <- sequences[as.data.frame(idx_match)[,1]]
  mis_idx_percent <- (init_len/length(sequences) - 1)*100

  # Grab barcode sequences
  # Trim 5' sequences
  sequences <- narrow(sequences, 1, 90)
  # Match pattern without accounting for indels
  match <- vmatchPattern(pat, sequences, length(strsplit(pat, 'N')[[1]]) + bb_mis)
  # Subset
  sequences_sub <- sequences[as.data.frame(match)[,1]]
  left_out <- sequences[-as.data.frame(match)[,1]]
  # Extract barcode sequences from sequences subset
  match <- vmatchPattern(pat, sequences_sub, length(strsplit(pat, 'N')[[1]]) + bb_mis)
  barcodes <- sequences_sub[match]
  if (indels) {
    # Remove sequences that will be matched by substitutions instead of indels (the total edit distance in vmatchPattern2 depends both on substitution and indels)
    substitution <- vmatchPattern(pat, left_out, max.mismatch = length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, min.mismatch = length(strsplit(pat, 'N')[[1]]) + bb_mis + indels) # length(strsplit(pat, 'N')[[1]]) returns the count of N in the barcode sequence (32)
    left_out <- left_out[-as.data.frame(substitution)[,1]]
    # Match pattern with remaining sequences, now accounting for indels
    match <- vmatchPattern2(pat, left_out, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T)
    # Subset
    indels_sub <- left_out[as.data.frame(match)[,1]]
    # Extract barcode sequences with indels
    match <- vmatchPattern2(pat, indels_sub, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T)
    barcodes <- c(barcodes, indels_sub[match])
  }

  # Discard sequences with 'N' calls
  barcodes <- barcodes[!grepl('N', barcodes)]
  
  # Generate count summary
  barcodes <- as.data.frame(barcodes)
  names(barcodes) <- 'seq'
  barcodes_summary <- count(barcodes, seq, sort=T)
  
  # Write summary file
  write.table(barcodes_summary, paste0(name, '_barcode_summary.txt'), sep='\t', row.names = F)
  write(mis_idx_percent, paste0(name, '_mis_idx.txt')) # This reports the percentage of sequences discarded by the index filtering
  
  # 1.b. Pool similar barcodes together (according to Hamming distance)
  #####################################################################
  
  message(paste('Pooling sample: ', name))
  barcodes_summary$seq <- as.character(barcodes_summary$seq)
  
  step <- 1
  while(step < (nrow(barcodes_summary))) {
    d <- stringdist(barcodes_summary$seq[step], barcodes_summary$seq[(step+1):length(barcodes_summary$seq)]) # calculate distance from top clone compared to the rest
    pool <- barcodes_summary[(step+1):length(barcodes_summary$seq),][d <= threshold,] # pool barcodes 
    barcodes_summary$n[step] <- barcodes_summary$n[step] + sum(pool$n) # update counts
    to_remove <- c((step+1):length(barcodes_summary$seq))[d <= threshold] # remove pooled barcodes from this step
    if(length(to_remove > 0)) {
      barcodes_summary <- barcodes_summary[-to_remove,]
    }
    step <- step + 1 # go to next line
  }
  
  # Sort summary table
  barcodes_summary <- arrange(barcodes_summary, desc(n))
  
  # Write summary file for pooled barcodes
  write.table(barcodes_summary, paste0(name, '_barcode_summary_pooled.txt'), sep='\t', row.names = F)
  message(paste('Finished sample: ', name))
}


# 2. Merge data
###############

message('Merging data')

# Load data summaries
summaries <- list.files()[grep('summary_pooled.txt', list.files())]
data.list <- list()
for (i in 1:length(summaries)) {
  data.list[[i]] <- read.table(summaries[i], sep = '\t', header = T)
}

# Merge data
merged <- NULL
current <- data.list[1]
for (i in 2:length(summaries)){
  merged <- merge(current, data.list[i], by = 'seq', all = T)
  current <- merged
}
colnames(merged) <- c('seq', as.character(strsplit(list.files()[grep('summary_pooled.txt', list.files())], '_barcode_summary_pooled.txt')))

# Replace NA's with 0
merged[is.na(merged)] <- 0

# Export merged table for counts
write.table(merged, 'merged_summary_pooled_count.txt', sep='\t', row.names = F)

# Convert counts to frequencies
for (i in 2:(length(summaries)+1)) {
  merged[,i] <- 100*merged[,i]/sum(merged[,i])
}

# Add the sum of frequencies at the last column of the merged data frame
sum.freq <- NULL
for (i in 1:nrow(merged)) {
  sum.freq <- c(sum.freq, sum(merged[i,2:(length(summaries)+1)]))
}
merged$sum.freq <- sum.freq
merged$seq <- as.character(merged$seq)

# Export merged table for frequencies
write.table(merged, 'merged_summary_pooled_freq.txt', sep='\t', row.names = F)

# export session info for reproducibility
writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
message('Done')