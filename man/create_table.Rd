% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_table.R
\name{create_table}
\alias{create_table}
\title{create a labeled coordinate table}
\usage{
create_table(im, threshold = 0, local_maxima_dist = 10,
  statname = "zstat", allow_duplicate_labels = TRUE, ...)
}
\arguments{
\item{im}{a\code{NeuroVol} instance}

\item{local_maxima_dist}{the distance threshold used to define local_maxima}

\item{statname}{the name of the statistic associated with the image values}

\item{allow_duplicate_labels}{whether to allow the same atlas label to appear more than once in the table. 
If \code{FALSE} and the same label appears twice, then only the coordinate with the highest vlaue will be retained.}

\item{...}{extra args to pass to \code{conn_comp}}

\item{a}{cluster threshold value to define the mask.}
}
\description{
create a labeled coordinate table
}
