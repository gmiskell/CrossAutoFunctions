#' A rolling Kolmogorov-Smirnov (KS) two-sample test function.
#'
#' This function gives the rolling KS test for an observation, relative to a proxy from other observations.
#' @param obs The time-series under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param proxy The comparison column.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @param window.length This defines the length of the sampled window, defined by row number.
#' @export
#' @examples
#' ksFUN()


ksFUN <- function(x, obs, date, proxy, reflective = TRUE, date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00', theta = NA, tau = NA, window.length = 72){
  
  # install and load required packages
    list.of.packages <- c("raster", "ggplot2", "plyr", "data.table", "stringr", "dplyr")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages) 
    library(raster); library(ggplot2); library(plyr); library(data.table); library(stringr); library(dplyr)
    
  # clause on type of analysis to be run
    if (reflective == TRUE) {
      
      date.start = date.start 
      date.end = date.end
    } else {
      
      date.start = str_c(Sys.Date()-7, ' 00:00:00')
      date.end = str_c(Sys.Date(), ' 00:00:00')
    }
    
    date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end)
  
  # define selected variables
  # use `data.table` package to deal with large datasets
    x <- as.data.table(x)
    setDT(x)
    x[, date := ymd_hms(date)]
    x <- x[date >= date.start & date <= date.end]
    x$proxy <- x[,..proxy]
    x$obs <- x[,..obs]
    
  # select the columns for day, groups, daily variance & n of sites
    x <- x[, list(date, obs, proxy)]
    
  # remove repeating and empty rows
    x <- unique.data.frame(x)
    x <- na.omit(x)

  # generate rolling KS test results using function by site
    ks.x <- rollingKStest(x, obs = 'obs', proxy = 'proxy', window = window.length)
  
  # use theta and tau thresholds if present
    if((!is.na(theta))) {
  
      if((!is.na(tau))) {

      # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
      # this looks at running means and is for tau, the length of time for alarms
      
        ks.x <- ks.x[, warning := ifelse(p.value < theta, 1, 0), by = list(group)][, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)] 
   
      } else {
    
        ks.x$alaram <- NA
      } 
    }  else {
  
      ks.x$warning <- NA
      ks.x$alarm <- NA
    }

  # gather variables of interest and return
    ks.x <- ks.x[, list(date, obs, proxy, test = 'ks test', warning, alarm)]

	ks.x
	
}
  