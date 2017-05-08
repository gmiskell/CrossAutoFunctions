#' A function for comparing the daily minima and maxima.
#'
#' This function finds and compares the daily 'minima' (5th) and 'maxima' (95th) percentiles across groups using modified z-scores (Mean Absolute Deviation - MAD).
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param group This is the first grouping value of the data, with the default being 'site'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to 2 (similar to 2 standard deviations).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' export
#' examples
#' percentFUN()

percentFUN <- function(x, date, obs, group, reflective = TRUE, date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00', theta = 2, tau = NA){
  
  # install and load required packages
    list.of.packages <- c("stats", "stringr", "reshape2", "ggplot2", "data.table", "dplyr")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    library(stats); library(stringr); library(reshape2); library(ggplot2); library(data.table); library(dplyr)
  
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
	  x <- data.table(x)
	  x$date <- setDT(x)[, ..date]
	  x[, date := ymd_hms(date)]
	  x <- x[date %between% c(date.start, date.end)]
	  x$group <- x[,..group]
	  x$obs <- x[,..obs]
  
  # add day column 
	  x <- x[, day := str_sub(date, end = 10)]
	
  # determine the 5th and 95th percentiles
    min.max <- x[, pm.5 := quantile(obs, probs = 0.05, na.rm = T), by = list(day, group)][, pm.95 := quantile(obs, probs = 0.95, na.rm = T), by = list(day, group)][, list(day, test = 'percentile', group, pm.5, pm.95)]

  # remove repeating and empty rows
	  min.max <- unique(min.max)
	  min.max <- na.omit(min.max)

  # derive daily MAD scores
	  min.max <- min.max[, min := median(pm.5, na.rm = T), by = day][, max := median(pm.95, na.rm = T), by = day][, min.sd := mad(pm.5, na.rm = T), by = day][, max.sd := mad(pm.95, na.rm = T), by = day]
  
  # determine MAD scores
	  min.max.test <- min.max[, min.test := (pm.5 - min)/(1.4826 * min.sd), by = list(day, group)][, max.test := (pm.95 - max)/(1.4826 * max.sd), by = list(day, group)]
	
  # use theta and tau thresholds
	  if((!is.na(theta))) {
  
		  if((!is.na(tau))) {
  
  # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
  # this looks at running means and is for tau, the length of time for alarms
        min.max.test <- min.max.test[, min.warning := ifelse(min.test > theta | min.test < -theta, 1, 0), by = list(group)][, min.alarm := movingFun(min.warning, n = tau, type = 'to', fun = mean, na.rm = T)][, max.warning := ifelse(max.test > theta | max.test < -theta, 1, 0), by = list(group)][, max.alarm := movingFun(max.warning, n = tau, type = 'to', fun = mean, na.rm = T)] 
		  } else {
	
    	  min.max.test$min.alarm <- NA; min.max.test$max.alarm <- NA
		  } 
	  } else {
	
      min.max.test$warning$min.warning <- NA; min.max.test$max.warning <- NA
	    min.max.test$alarm$min.alarm <- NA; min.max.test$max.alarm <- NA
	  }
  
  # gather variables of interest and return
	  min.max.test <- min.max.test[, list(date = day, group, test = 'percent', min.warning, min.alarm, max.warning, max.alarm)]
	     
	return(min.max.test)
}