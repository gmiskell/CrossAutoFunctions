#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD'.
#' @param group1 This is the first grouping value of the data, with the default being 'site'.
#' @param group2 This is the second grouping value of the data, which allows overlain networks, with the default being 'pol'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param dest This is the destination folder where the function output will be saved. Defaults to current working directory.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' require(lubridate)
#' data(lakers)
#' rvFUN(lakers, obs = 'x', group1 = 'period', group2 = 'team', date.start = '2008-10-28', date.end = '2009-04-14')

rvFUN <- function(x, obs, group1, group2, reflective = TRUE, plot = TRUE, dest = getwd(), date.start = '2016-07-01', date.end = '2016-07-10', theta = NA, tau = NA) {
  
  # install and load required packages
  list.of.packages <- c("stats", "stringr", "raster", "ggplot2", "data.table", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(stats); library(stringr); library(raster); library(ggplot2); library(data.table); library(dplyr)
  
  # clause on type of analysis to be run, define dates
	if (reflective == TRUE) {
	
		date.start = date.start 
		date.end = date.end
	} else {
	
		date.start = str_c(Sys.Date()-7)
		date.end = str_c(Sys.Date())
	}
	
  # define selected variables
  # use data.table package to deal with large datasets
	setDT(x)
	x[, date := ymd(date)]
	x <- x[date >= date.start & date <= date.end]
	x$group1 <- x[,..group1]
	x$group2 <- x[,..group2]
	x$obs <- x[,..obs]
	
  # add day column 
	x <- x[, day := str_sub(date, end = 10)]
  
  # add daily variance and daily number of group1 across day and group2
	x <- x[, day.var := var(obs, na.rm = T), by = list(group1, group2, day)]
	x <- x[, n.grp := length(unique(group1)), by = list(day, group2)]
  
  # select the columns for day, groups, daily variance & n of sites
	x <- x[, list(day, test = 'rank variance', group1, group2, day.var, n.grp)]
  
  # remove repeating and empty rows
	x <- unique.data.frame(x)
	x <- na.omit(x)
  
  # create ecdf column (this will create warnings)
	x <- suppressWarnings(x[, result := ecdf(day.var)(unique(day.var)), by = list(day, group2)])
  
  # use theta and tau thresholds if present
	if((!is.na(theta))) {
  
		if((!is.na(tau))) {
  
			if (theta < 0.5) {
    
  # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
  # this looks at running means and is for tau, the length of time for alarms
	
			x <- x[, change.thresh := ifelse(result < theta, 1, 0), by = list(group1, group2)][, alarm := movingFun(change.thresh, n = tau, type = 'to', fun = mean, na.rm = T)] 
			} else {   

  # this is similar to when theta < 0.5
	        x <- x[, change.thresh := ifelse(result > theta, 1, 0), by = list(group1, group2)][, alarm := movingFun(change.thresh, n = tau, type = 'to', fun = mean, na.rm = T)]
			} 
		} else {
		
		tau <- NA
		} 
	}  else {
	
	theta <- NA
	tau <- NA
	}
  
  # create csv file to be saved to dest folder (default is directory)
  csv <- write.csv(x, file = str_c(dest, '/rankedVariance', start.date, '.csv'))
  
  # create plot and save
  if(plot == TRUE) {
  
  # create plot of the ranked variance over the last week
	plot.name <- str_c(dest, '/rankedVariance', start.date, '.pdf')
	pdf(plot.name, width = 6, height = 10)
	ts <- suppressWarnings(ggplot(x, aes(x = as.POSIXct(day), y = result, colour = group2)) +
				theme_bw() +
				geom_point() +
				facet_wrap(~ group1) +
				xlab("") + ylab('ECDF'))
  
		if(!is.na(theta)) {ts <- ts + geom_abline(intercept = theta, slope = 0)}
	dev.off()
	return(c(ts, csv))
  } else { 
  
  # returns x as saved csv file
  csv
  }
}