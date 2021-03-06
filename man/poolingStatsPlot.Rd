% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poolingStatsPlot.R
\name{poolingStatsPlot}
\alias{poolingStatsPlot}
\title{A function to export QC plots after the pooling step.}
\usage{
poolingStatsPlot(dir = getwd())
}
\arguments{
\item{dir}{The directory containing the statistics from the pooling step.}
}
\description{
A function to export QC plots after the pooling step.
}
\details{
This function shows the frequency of the pooled barcodes according to the distance to the reference barcode.
}
\examples{
poolingStatsPlot(dir = getwd())
}
\keyword{BC32}
\keyword{pooling}
\keyword{processing}
\keyword{stats}
