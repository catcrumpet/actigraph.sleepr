% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_epochs.R
\name{collapse_epochs}
\alias{collapse_epochs}
\title{Re-integrate epochs}
\usage{
collapse_epochs(agdb, epoch_len_out, use_incomplete = TRUE)
}
\arguments{
\item{agdb}{A \code{tibble} (\code{tbl}) of activity data (at least) an \code{epochlength} attribute.}

\item{epoch_len_out}{Output (longer) epoch length in seconds, must be exact multiple of the input epoch length.}

\item{use_incomplete}{logical. Set to \code{TRUE} to follow ActiLife convention, which collapses all observed epochs even if they are incomplete.}
}
\value{
A \code{tibble} (\code{tbl}) of activity data collapsed into epochs of specified length.
}
\description{
Collapse post-filtered activity counts into larger epoch "buckets".
}
\details{
Output epochs start at the first timestamp of the input data.
}
\examples{
library("dplyr")
data("gtxplus1day")

gtxplus1day \%>\%
  collapse_epochs(30)

gtxplus1day \%>\%
  collapse_epochs(60)
}
\references{
ActiLife 6 User's Manual by the ActiGraph Software Department. 04/03/2012.
}
