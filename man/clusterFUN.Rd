% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterFUN.r
\name{clusterFUN}
\alias{clusterFUN}
\title{A daily cluster function.}
\usage{
clusterFUN(x, date, obs, group, by.day = TRUE, reflective = TRUE,
  theta = NA, tau = NA)
}
\arguments{
\item{date}{A date column (set up as `YYYY-MM-DD HH:MM:SS`).}

\item{obs}{The column under observation.}

\item{group}{A column which the data are partitioned on (either spatial or temporal).}

\item{by.day}{An option where analysis uses data from only that day ("cross", `TRUE`, default), or uses all selected data ("auto", `FALSE`). If selecting auto, the last five days will be used.}

\item{reflective}{Decides whether to use a selected time period of data (reflective, `TRUE`, default), or to use data from the latest day (live, `FALSE`).}

\item{theta}{Sets the test threshold on whether to flag the data. Defaults to `NA` (no flags given).}

\item{tau}{Sets the day threshold on whether to flag the data given consistent theta flags. Defaults to `NA` (no flags given).}
}
\description{
A daily cluster function.
}
\examples{
clusterFUN()
}
