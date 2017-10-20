#' A function to generate summaries from BC32 fastq files.
#'
#' @param pat The barcode backbone pattern for matching.
#' @param restriction The nucleotide sequence of the restriction site flanking the barcode at the 3' end.
#' @param sampname A data frame with 3 columns containing the name of each sample to
#' process, its file name and the expected multiplex index sequence.
#' @param base_q The minimum mean base quality in the first 90 nucleotides required to keep
#' a read for the processing. Defaults to 20.
#' @param idx_mis The number of mismatches allowed during index matching. Defaults to 1.
#' @param bb_mis The number of mismatches allowed during barcode matching. Defaults to 1.
#' @param indels The number of indels allowed during barcode matching. This strongly impacts
#' the speed of the function. Try to keep it low. Defaults to 1.
#' @details This function extracts the number of different barcode sequences found in a provided
#' list of fastq files and export a summary table, quality control plots, and a list of the most
#' common sequences that failed index or barcode matching.
#' @return Returns data in a nested list of samples containing a list of sequences not matching the index,
#' QC data and sequences not matching the barcode index.
#' @keywords BC32 processing QC
#' @export
#' @examples
#' bc_data <- generateSummaries(pat = "CTAGCCAGTT",
#'                              restriction = "CTCGAG"
#'                              sampname = sampname)

generateSummaries <- function(pat,
                              restriction,
                              sampname,
                              base_q = 20,
                              idx_mis = 1,
                              bb_mis = 1,
                              indels = 1) {

  to_return <- list() # returned object

  # Loop through samples
  for (i in 1:length(sampname$file)) {

    # Read fastq file with ShortRead package
    sequences <- readFastq(sampname$file[i])
    init_len <- length(sequences) # collect for QC_stats
    name <- sampname$sample[i]
    message(paste('Processing sample: ', name))

    # Filter reads with average quality below base_q
    good_qual <- (alphabetScore(quality(narrow(sequences, 1, 90)))/90) >= base_q
    sequences <- sequences[good_qual]@sread
    good_reads <- length(sequences)

    # Filter for sequences with specific multiplexing index (only keep sequences with the correct index, beware of 6 and 8-mer indexes!)
    # This steps avoids misassignment of samples at de-multiplexing
    message('index matching')
    idx <- paste0('GTCAC', sampname$index[i], 'ATCTC')
    # Index matching
    idx_match <- vmatchPattern(idx, sequences, idx_mis)
    # Subset sequences
    if (length(as.data.frame(idx_match)[,1]) > 0) {
      no_index_seq <- sequences[-as.data.frame(idx_match)[,1]] # keep sequences failing index matching
    } else {
      no_index_seq <- sequences
    }
    to_return[[i]] <- list(no_index_seq)
    sequences <- sequences[as.data.frame(idx_match)[,1]]
    correct_index <- length(sequences) # collect for QC_stats

    # Grab barcode sequences
    message('barcode matching')
    # Trim 3' end
    sequences <- narrow(sequences, 1, 90)
    # Match pattern without accounting for indels
    match <- vmatchPattern(pat, sequences, bb_mis, fixed = FALSE)
    # Subset
    sequences_sub <- sequences[as.data.frame(match)[,1]]
    match_no_indels <- length(sequences_sub) # collect for QC_stats
    left_out <- sequences[-as.data.frame(match)[,1]]
    # Extract barcode sequences from sequences subset
    match <- vmatchPattern(pat, sequences_sub, bb_mis, fixed = FALSE)
    barcodes <- sequences_sub[match]
    match_with_indels <- 'sample processed without indel matching'
    if (indels) {
      message('accounting for indels')
      pat_indels <- paste0(pat, restriction) # include restriction sequence in the matching pattern
      # Remove sequences that will be matched by additional substitutions instead of indels (the total edit distance in vmatchPattern2 depends both on substitution and indels)
      substitution <- vmatchPattern(pat_indels, left_out, bb_mis + indels, fixed = FALSE)
      if (length(as.data.frame(substitution)[,1]) > 0) {
        left_out <-  left_out[-as.data.frame(substitution)[,1]]
      }
      # Match pattern with remaining sequences, now accounting for indels
      match <- vmatchPattern2(pat_indels, left_out, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T) # length(strsplit(pat, 'N')[[1]]) returns the count of N in the barcode sequence (32)
      # Subset and remove sequences with too long edit distance
      to_include <- as.data.frame(match)[as.data.frame(match)[,5] >= nchar(pat_indels) - indels,][,1]
      indels_sub <- left_out[to_include]
      if (length(to_include) > 0) {
        left_out <- left_out[-to_include]
      }
      # Extract barcode sequences with indels
      match <- vmatchPattern2(pat_indels, indels_sub, length(strsplit(pat, 'N')[[1]]) + bb_mis + indels, with.indels = T)
      indels_sub <- indels_sub[match]
      indels_sub <- trimLRPatterns(Rpattern = restriction, subject = indels_sub, max.Rmismatch = indels, with.Rindels = T) # trim restriction sequence pattern
      match_with_indels <- length(indels_sub) # collect for QC_stats
      barcodes <- c(barcodes, indels_sub)
    }

    # Discard sequences with 'N' calls
    barcodes <- barcodes[!grepl('N', barcodes)]

    # Generate count summary
    message('generating outputs')
    barcodes <- as.data.frame(barcodes)
    names(barcodes) <- 'seq'
    barcodes_summary <- dplyr::count(barcodes, seq, sort=T)

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
    to_return[[i]][[length(to_return[[i]]) + 1]] <- QC_stats
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

    # Export the top-10 reads not matching index or barcode to investigate causes (also print out alignment for help)

    no_index_seq <- as.data.frame(head(dplyr::count(as.data.frame(no_index_seq), x, sort=T), n = 10))
    aln <- NULL
    for (j in 1:nrow(no_index_seq)){
      aln <- c(aln, msa(c(idx, no_index_seq[j,1]), method = 'ClustalOmega', type = 'dna'))
    }
    write.table(c(capture.output(no_index_seq), capture.output(aln)), paste0(name, '_no_index.txt'), sep = '\t', row.names = F, col.names = F, quote = F)

    to_return[[i]][[length(to_return[[i]]) + 1]] <- as.data.frame(left_out)
    no_bc_match <- as.data.frame(head(dplyr::count(as.data.frame(left_out), x, sort=T), n = 10))
    aln <- NULL
    for (j in 1:nrow(no_bc_match)){
      aln <- c(aln, msa(c(pat, no_bc_match[j,1]), method = 'ClustalOmega', type = 'dna'))
    }
    write.table(c(capture.output(no_bc_match), capture.output(aln)), paste0(name, '_no_bc_match.txt'), sep = '\t', row.names = F, col.names = F, quote = F)
    names(to_return[[i]]) <- c('no_index', 'QC_stats', 'no_bc')
  }
  names(to_return) <- sampname$sample
  message(paste('Sample ', name, ' processed'))
  to_return
}
