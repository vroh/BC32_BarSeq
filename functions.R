### Helping functions for Barcode Processing script ###

# Below is an updated version of vmatchPattern that take indels into account
# Provided by Hervé Pagès on https://support.bioconductor.org/p/58350/
vmatchPattern2 <- function(pattern, subject,
                           max.mismatch=0, min.mismatch=0,
                           with.indels=FALSE, fixed=TRUE,
                           algorithm="auto")
{
  if (!is(subject, "XStringSet"))
    subject <- Biostrings:::XStringSet(NULL, subject)
  algo <- Biostrings:::normargAlgorithm(algorithm)
  if (Biostrings:::isCharacterAlgo(algo))
    stop("'subject' must be a single (non-empty) string ",
         "for this algorithm")
  pattern <- Biostrings:::normargPattern(pattern, subject)
  max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
  min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch,
                                                  max.mismatch)
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  algo <- Biostrings:::selectAlgo(algo, pattern,
                                  max.mismatch, min.mismatch,
                                  with.indels, fixed)
  C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
                  max.mismatch, min.mismatch,
                  with.indels, fixed, algo,
                  "MATCHES_AS_RANGES",
                  PACKAGE="Biostrings")
  unlisted_ans <- IRanges(start=unlist(C_ans[[1L]],
                                       use.names=FALSE),
                          width=unlist(C_ans[[2L]],
                                       use.names=FALSE))
  relist(unlisted_ans, C_ans[[1L]])
}