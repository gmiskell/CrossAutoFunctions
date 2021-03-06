% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ksFUN.R
\name{ksFUN}
\alias{ksFUN}
\title{A rolling Kolmogorov-Smirnov (KS) two-sample test function.}
\usage{
ksFUN(x, obs, date, proxy, reflective = TRUE, theta = NA, tau = NA,
  window.length = 72)
}
\arguments{
\item{obs}{The time-series under investigation.}

\item{date}{This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.}

\item{proxy}{The comparison column.}

\item{reflective}{This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.}

\item{theta}{This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).}

\item{tau}{This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).}

\item{window.length}{This defines the length of the sampled window, defined by row number.}

\item{date.start}{This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.}

\item{date.end}{This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.}
}
\description{
This function gives the rolling KS test for an observation, relative to a proxy from other observations.
}
\examples{
ksFUN()
}
