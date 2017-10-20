#' String searching function.
#'
#' A functions for finding all the occurrences (aka "matches" or "hits") of a given pattern (typically
#' short) in a (typically long) reference sequence or set of reference sequences (aka the subject).
#' This is an updated version of vmatchPattern that take indels into account.
#' Provided by Hervé Pagès on https://support.bioconductor.org/p/58350/
#' @param pattern The pattern string.
#' @param subject The String object for matching.
#' @param max.mismatch,min.mismatch The minimum and maximum number of mismatching letters allowed.
#' @param with.indels If TRUE then indels are allowed.
#' @param fixed If TRUE (the default), an IUPAC ambiguity code in the pattern can only match
#' the same code in the subject, and vice versa.
#' @param algorithm One of the following: "auto", "naive-exact", "naive-inexact", "boyer-moore",
#' "shift-or" or "indels".
#' @param ... Additional arguments for methods.
#' @details See vmatchPattern for details.
#' @keywords pattern matching
#' @export
#' @examples
#' x <- DNAString("AAGCGCGATATG")
#' m1 <- matchPattern("GCNNNAT", x)
#' m1

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
