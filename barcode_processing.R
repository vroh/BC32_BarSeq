##################################################
#                                                #
#          Barcode Processing v.0.4              #
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
library(ggplot2)
library(gridExtra)

# For indels, an updated version of Biostrings' vmatchPattern is required (described in functions.R)
source('functions.R')

# Define variables
pat <- 'CTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNN' # matching pattern (BC32)
base_q <- 0 # minimum average base quality score in first 90 nucleotides
bb_mis <- 1 # number of mismatch allowed in barcode backbone sequence
indels <- 0 # total edit distance resulting from indels that are tolerated for barcode matching (try to keep this low)
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
  sequences <- readFastq(seqfiles[i])
  init_len <- length(sequences) # collect for QC_stats
  name <- sampname$sample[i]
  message(paste('Processing sample: ', name))
  
  # Filter reads with average quality below base_q
  good_qual <- (alphabetScore(quality(narrow(sequences, 1, 90)))/90) >= base_q
  sequences <- sequences[good_qual]@sread
  good_reads <- length(sequences)
  
  # Filter for sequences with specific multiplexing index (only keep sequences with the correct index, beware of 6 and 8-mer indexes!)
  # This steps avoids misassignment of samples at de-multiplexing 
  idx <- paste0('GTCAC', sampname$index[i], 'ATCTC')
  # Index matching (strict!)
  idx_match <- vmatchPattern(idx, sequences, 0)
  # Subset sequences
  sequences <- sequences[as.data.frame(idx_match)[,1]]
  correct_index <- length(sequences) # collect for QC_stats

  # Grab barcode sequences
  # Trim 5' sequences
  sequences <- narrow(sequences, 1, 90)
  # Match pattern without accounting for indels
  match <- vmatchPattern(pat, sequences, length(strsplit(pat, 'N')[[1]]) + bb_mis)
  # Subset
  sequences_sub <- sequences[as.data.frame(match)[,1]]
  match_no_indels <- length(sequences_sub) # collect for QC_stats
  left_out <- sequences[-as.data.frame(match)[,1]]
  # Extract barcode sequences from sequences subset
  match <- vmatchPattern(pat, sequences_sub, length(strsplit(pat, 'N')[[1]]) + bb_mis)
  barcodes <- sequences_sub[match]
  match_with_indels <- 'sample processed without indel matching'
  if (indels) {
    pat_indels <- paste0(pat, 'CTCGAG') # include 'CTCGAG' in the matching pattern
    # Remove sequences that will be matched by additional substitutions instead of indels (the total edit distance in vmatchPattern2 depends both on substitution and indels)
    substitution <- vmatchPattern(pat_indels, left_out, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels) # length(strsplit(pat, 'N')[[1]]) returns the count of N in the barcode sequence (32)
    if (length(as.data.frame(substitution)[,1]) > 0) {
      left_out <-  left_out[-as.data.frame(substitution)[,1]]
    }
    # Match pattern with remaining sequences, now accounting for indels
    match <- vmatchPattern2(pat_indels, left_out, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T)
    # Subset and remove sequences with too long edit distance
    indels_sub <- left_out[as.data.frame(match)[as.data.frame(match)[,5] >= nchar(pat_indels) - indels,][,1]]
    # Extract barcode sequences with indels
    match <- vmatchPattern2(pat_indels, indels_sub, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T)
    indels_sub <- indels_sub[match]
    indels_sub <- trimLRPatterns(Rpattern = 'CTCGAG', subject = indels_sub, max.Rmismatch = indels, with.Rindels = T) # trim CTCGAG pattern
    match_with_indels <- length(indels_sub) # collect for QC_stats
    barcodes <- c(barcodes, indels_sub)
  }

  # Discard sequences with 'N' calls
  barcodes <- barcodes[!grepl('N', barcodes)]
  
  # Generate count summary
  barcodes <- as.data.frame(barcodes)
  names(barcodes) <- 'seq'
  barcodes_summary <- count(barcodes, seq, sort=T)
  
  # Write summary file and QC_stats report
  write.table(barcodes_summary, paste0(name, '_barcode_summary.txt'), sep='\t', row.names = F)
  QC_stats <- rbind(c('Number of FASTQ sequences', init_len),
                    c('Number of reads above quality threshold', good_reads),
                    c('% good quality reads', 100*(good_reads/init_len)),
                    c('Sequences matching multiplex index', correct_index),
                    c('% correct index', 100*(correct_index/good_reads)),
                    c('Barcodes match (without indels)', match_no_indels),
                    c('% without indels', 100*(match_no_indels/correct_index)),
                    c('Barcodes matching (with indels)', match_with_indels),
                    c('% with and without indels', 100*(match_no_indels + ifelse(indels, match_with_indels, 0))/correct_index))
  write.table(QC_stats, paste0(name, '_QC_stats.xls'), sep = '\t', row.names = F, col.names = F)
  
  # Export QC plots
  # prepare data frame
  QC_plot <- as.data.frame(QC_stats[c(1,2,4),])
  QC_plot$V2 <- as.numeric(as.character(QC_plot$V2))
  # export plots as pdf
  pdf(paste0(name, '_QC_plots.pdf'), 10)
  # plotA: sequences after quality filtering and index matching
  plotA <- ggplot(QC_plot, aes(V1, V2)) +
    geom_bar(stat = 'identity') +
    ylim(0, max(QC_plot$V2)) +
    xlab("") +
    ylab("Count") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_label(aes(x = 1, y = 0, label = paste0('n = ', QC_stats[1,2]))) +
    geom_label(aes(x = 2, y = 0, label = paste0('% = ', round(as.numeric(QC_stats[3,2]), 2)))) +
    geom_label(aes(x = 3, y = 0, label = paste0('% = ', round(as.numeric(QC_stats[5,2]), 2)))) +
    geom_label(aes(x = 1, y = max(QC_plot$V2), label = 'FASTQ sequences')) +
    geom_label(aes(x = 2, y = max(QC_plot$V2), label = 'Good quality')) +
    geom_label(aes(x = 3, y = max(QC_plot$V2), label = 'Correct index'))
  # plotB: sequences after barcode matching
  QC_plot_indels <- as.data.frame(cbind(c('substitutions', 'indels'),
                                        c(100*(as.numeric(as.character(QC_stats[6,2]))/as.numeric(as.character(QC_stats[4,2]))),
                                          100*(as.numeric(as.character(QC_stats[8,2]))/as.numeric(as.character(QC_stats[4,2]))))))
  QC_plot_indels$V2 <- as.numeric(as.character(QC_plot_indels$V2))
  plotB <- ggplot(QC_plot_indels, aes(1, V2)) +
    geom_bar(stat = 'identity', aes(fill = V1)) +
    ylim(0, 100) +
    xlab("") +
    ylab("% sequences") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_fill_discrete(guide_legend(title = 'matching group')) +
    geom_label(aes(x = 1, y = 0, label = paste0('% = ', round(as.numeric(QC_stats[9,2]), 2)))) +
    geom_label(aes(x = 1, y = 100, label = 'Barcode matching'))
  grid.arrange(plotA, plotB, nrow = 1)
  dev.off()
  
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
