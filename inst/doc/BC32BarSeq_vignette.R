## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  library(BC32BarSeq)

## ---- eval=FALSE---------------------------------------------------------
#  initBC()

## ---- eval=FALSE---------------------------------------------------------
#  pat <- 'CTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNNCTTNNCGANNCTANNCTTNNGGANNCTANNCAGNN' # matching pattern (BC32)
#  restriction <- 'CTCGAG' # The restriction sequence flanking the barcode pattern at the 3'end
#  idx_mis <- 1 # number of mismatch allowed in index
#  base_q <- 20 # minimum average base quality score in first 90 nucleotides
#  bb_mis <- 1 # number of mismatch allowed in barcode backbone sequence
#  indels <- 1 # total edit distance resulting from indels that are tolerated for barcode matching (try to keep this low)
#  dl <- 12 # threshold for Damerau-Levenshtein distance in step pooling similar sequences
#  min.count <- 10 # in the final merged table, replace barcodes with occurence below min.count (considered as noise) with 0

## ---- eval=FALSE---------------------------------------------------------
#  sampname <- read.delim('sampname.txt', # sampname provides a list of sample names, matching sequencing files and multiplexing index
#                         header = F,
#                         stringsAsFactor = F)
#  colnames(sampname) <- c('sample', 'file', 'index')

## ---- eval=FALSE---------------------------------------------------------
#  bc_data <- generateSummaries(pat = pat,
#                               restriction = restriction,
#                               sampname = sampname,
#                               base_q = base_q,
#                               idx_mis = idx_mis,
#                               bb_mis = bb_mis,
#                               indels = indels)
#  
#  save(bc_data, 'bc_data.Rdata) # Keep this object, it is required for downstream processing!

## ---- eval=FALSE---------------------------------------------------------
#  invited_idx <- read.delim('invited_idx.txt', header = F, stringsAsFactor = F) # format: column 1 -> sequences, column 2 -> sample name
#  invited_bc <- read.delim('invited_bc.txt', header = F, stringsAsFactor = F) # format: column 1 -> sequences
#  
#  addInvitedSeq(pat = pat,
#                bc_data = bc_data,
#                dir = getwd(),
#                invited_idx = invited_idx,
#                invited_bc = invited_bc)

## ---- eval=FALSE---------------------------------------------------------
#  poolBC(sampname = sampname,
#          dir = getwd(),
#          dl = dl)

## ---- eval=FALSE---------------------------------------------------------
#  counts <- mergeSummaries(sampname = sampname,
#                           dir = getwd(),
#                           min.count = min.count)
#  
#  freq <- countsToFreq(counts = counts)
#  
#  cpm <- countsToCPM(counts = counts)
#  
#  write.table(counts, 'merged_summary_pooled_count.txt', sep='\t', row.names = F)
#  write.table(freq, 'merged_summary_pooled_freq.txt', sep='\t', row.names = F)
#  write.table(cpm, 'merged_summary_pooled_cpm.txt', sep='\t', row.names = F)
#  
#  # export session info for reproducibility
#  writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

## ---- eval=FALSE---------------------------------------------------------
#  exploreBC(getwd())

