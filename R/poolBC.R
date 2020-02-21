#' A function to pool similar barcodes according to Damerau-Levenshtein distance.
#'
#' @param sampname A data frame with 3 columns containing the name of each sample to
#' process, its file name and the expected multiplex index sequence.
#' @param dir The directory containing the summaries.
#' @param dl The maximum Damerau-Levenshtein distance to consider barcodes as identical and pool them together. Defaults to 12.
#' @details This function updates the count summaries by pooling together similar barcodes (according to the given Damerau-Levenshtein distance).
#' @keywords BC32 processing update pool Damerau-Levenshtein
#' @export
#' @examples
#' poolBC(sampname, getwd(), 12)

poolBC <- function(sampname,
                    dir, # directory containing the summaries
                    dl = 12) {

  for (i in 1:length(bc_data)) {
    name <- sampname$sample[i]
    message(paste0('pooling sample: '), name)
    barcodes_summary <- read.delim(paste0(dir, '/', name, '_barcode_summary.txt'), stringsAsFactors = F)
    sample_depth <- sum(barcodes_summary$n)
    pooling_stats <- NULL

    step <- 1
    while(step < (nrow(barcodes_summary))) {
      stat_current <- NULL # create object that will store stats
      d <- stringdist(barcodes_summary$seq[step], barcodes_summary$seq[(step+1):length(barcodes_summary$seq)]) # calculate distance from top clone compared to the rest
      matching_vector <- d <= dl # keep track of barcodes to pool
      pool <- barcodes_summary[(step+1):length(barcodes_summary$seq),][matching_vector,] # pool barcodes

      # record stats
      stat_current$alt <- pool$seq # record alternative sequences
      stat_current$ref <- barcodes_summary$seq[step] # get name of current reference sequence
      stat_current$distance <- d[matching_vector] # record distances
      stat_current$rank <- (1:length(matching_vector))[matching_vector] # record ranks of matching barcodes
      stat_current$alt_counts <- pool$n
      stat_current$alt_freq <- pool$n / sample_depth
      stat_current$ref_counts <- barcodes_summary$n[step]
      stat_current$ref_freq <- barcodes_summary$n[step] / sample_depth

      # continue pooling
      barcodes_summary$n[step] <- barcodes_summary$n[step] + sum(pool$n) # update counts
      to_remove <- c((step+1):length(barcodes_summary$seq))[d <= dl] # remove pooled barcodes from this step
      if(length(to_remove > 0)) {
        barcodes_summary <- barcodes_summary[-to_remove,]
      }
      if(nrow(pool) > 0) {
        pooling_stats <- rbind(pooling_stats, as.data.frame(stat_current))
      }
      step <- step + 1 # go to next line
    }

    # Sort summary table
    barcodes_summary <- arrange(barcodes_summary, desc(n))

    # Write summary file and stat file for pooled barcodes
    write.table(barcodes_summary, paste0(dir, '/', name, '_barcode_summary_pooled.txt'), sep='\t', row.names = F)
    write.table(pooling_stats, paste0(dir, '/', name, '_pooling_stats.txt'), sep='\t', row.names = F)

  }
}
