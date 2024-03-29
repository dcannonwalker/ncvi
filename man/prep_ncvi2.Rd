% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_ncvi2.R
\name{prep_ncvi2}
\alias{prep_ncvi2}
\title{Prepare list of objects required for \code{fit_ncvi()}. This is an
updated version of \code{prep_ncvi()}, which will replace it once
I've done some debugging.}
\usage{
prep_ncvi2(counts, S, X, Z)
}
\arguments{
\item{counts}{Data frame of counts with first column the
gene ID}

\item{S}{Vector of size factors with \code{length(S) = ncol(counts) - 1}}

\item{X}{Fixed effects design matrix (expects preparation in column 2
and treatment in column 3)}

\item{Z}{Random effects design matrix}
}
\description{
Prepare list of objects required for \code{fit_ncvi()}. This is an
updated version of \code{prep_ncvi()}, which will replace it once
I've done some debugging.
}
