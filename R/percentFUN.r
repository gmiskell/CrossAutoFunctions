#' A function for comparing the minima and maxima
#'
#' This function finds and compares the daily 'minima' (5th) and 'maxima' (95th) percentiles across groups using modified z-scores (Mean Absolute Deviation - MAD).
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD'.
#' @param group1 This is the first grouping value of the data, with the default being 'site'.
#' @param group2 This is the second grouping value of the data, which allows overlain networks, with the default being 'pol'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param dest This is the destination folder where the function output will be saved. Defaults to current working directory.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to 2 (similar to 2 standard deviations).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' export
#' examples
#' require(lubridate)
#' data(lakers)
#' percentFUN(lakers, obs = 'x', group1 = 'period', group2 = 'team', date.start = '2008-10-28', date.end = '2009-04-14')

percentFUN <- function(x, obs, group1, group2, reflective = TRUE, plot = TRUE, dest = getwd(), date.start = '2016-07-01', date.end = '2016-07-10', theta = 2, tau = NA){
  
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
	
		date.start = str_c(Sys.Date()-7)
		date.end = str_c(Sys.Date())
	}
  
  # define selected variables
  # use package to deal with large datasets
	setDT(x)
	x[, date := ymd(date)]
	x <- x[date >= date.start & date <= date.end]
	x$group1 <- x[,..group1]
	x$group2 <- x[,..group2]
	x$obs <- x[,..obs]
  
  # add day column 
	x <- x[, day := str_sub(date, end = 10)]
	
  # determine the 5th and 95th percentiles
    min.max <- x[, pm.5 := quantile(obs, probs = 0.05, na.rm = T), by = list(day, group1)][, pm.95 := quantile(obs, probs = 0.95, na.rm = T), by = list(day, group1)][, list(day, test = 'percentile', group1, group2, pm.5, pm.95)]

  # remove repeating and empty rows
	min.max <- unique.data.frame(min.max)
	min.max <- na.omit(min.max)

  # derive daily MAD scores
	min.max <- min.max[, min := median(pm.5, na.rm = T), by = day][, max := median(pm.95, na.rm = T), by = day][, min.sd := mad(pm.5, na.rm = T), by = day][, max.sd := mad(pm.95, na.rm = T), by = day]
  
  # determine MAD scores
	min.max.test <- min.max[, min.test := (pm.5 - min)/(1.4826 * min.sd), by = list(day, group1, group2)][, max.test := (pm.95 - max)/(1.4826 * max.sd), by = list(day, group1, group2)]
	
  # use theta and tau thresholds
	if((!is.na(theta))) {
  
		if((!is.na(tau))) {
  
  # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
  # this looks at running means and is for tau, the length of time for alarms
        min.max.test <- min.max.test[, min.thresh := ifelse(min.test > theta | min.test < -theta, 1, 0), by = list(group1, group2)][, min.alarm := movingFun(min.thresh, n = tau, type = 'to', fun = mean, na.rm = T)][, max.thresh := ifelse(max.test > theta | max.test < -theta, 1, 0), by = list(group1, group2)][, max.alarm := movingFun(max.thresh, n = tau, type = 'to', fun = mean, na.rm = T)] 
		} else {
	
    	tau <- NA
		} 
	}  else {
	
	theta <- NA
	tau <- NA
	}
  
  # write csv to dest folder (default is directory)
  csv <- write.csv(min.max.test, file = str_c(dest, '/minmaxPercent', start.date, '.csv'))
  
  if(plot == TRUE) {
  
  # create plot of the ranked variance over the last week
	plot.name <- str_c(dest, '/minmaxPercent', start.date, '.pdf')
    pdf(plot.name, width = 6, height = 10)
    min.ts <- suppressWarnings(ggplot(min.max.test, aes(x = as.POSIXct(day), y = min.test, colour = group2)) +
				theme_bw() +
				geom_line() +
				facet_wrap(~ group1) +
				xlab("") + ylab('daily minima'))
	max.ts <- suppressWarnings(ggplot(min.max.test, aes(x = as.POSIXct(day), y = max.test, colour = group2)) +
				theme_bw() +
				geom_line() +
				facet_wrap(~ group1) +
				xlab("") + ylab('daily maxima'))
		if(!is.na(theta)) { min.ts <- min.ts + geom_abline(intercept = c(-theta, theta), slope = 0); max.ts <- max.ts + geom_abline(intercept = c(-theta, theta), slope = 0)}
	dev.off()
    return(c(ts, min.ts, max.ts))
  } else { 
  
  # returns x as saved csv file
  csv
  }  
}