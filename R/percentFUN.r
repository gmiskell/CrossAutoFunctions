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
	
    date.start = Sys.time()-60*60*24*7
    date.end = Sys.time()
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
  min.test <- x[, quart := quantile(obs, probs = 0.05, na.rm = T), by = list(day, group)][, list(day, test = 'min_percentile', group, quart)]
  max.test <- x[, quart := quantile(obs, probs = 0.95, na.rm = T), by = list(day, group)][, list(day, test =  'max_percentile', group, quart)]
   min.max <- rbind(min.test, max.test)

  # remove repeating rows
  min.max <- unique(min.max)

  # derive daily MAD scores
  min.max <- min.max[, median := median(quart, na.rm = T), by = list(day, test)][, sd := mad(quart, na.rm = T), by = list(day, test)]
  
  # determine MAD score effect sizes
  min.max.test <- min.max[, sigma := (quart - median)/(1.4826 * sd)]
	
  # use theta and tau thresholds
  if((!is.na(theta))) {
  
	if((!is.na(tau))) {
  
    # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
    # this looks at running means and is for tau, the length of time for alarms
		min.max.test <- min.max.test[, warning := ifelse(abs(sigma) > theta, 1, 0)][, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T), by = list(group, test)]
	} else {
	
		min.max.test$alarm <- NA
	} 
  } else {
	
    min.max.test$warning <- NA
	min.max.test$alarm <- NA
  }
  
  # gather variables of interest and return
  min.max.test <- min.max.test[, list(date = as.POSIXct(day), site = group, test, statistic = sigma, warning, alarm)]
  return(min.max.test)

}
