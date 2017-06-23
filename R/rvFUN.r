#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param group This is the first grouping value of the data, with the default being 'site'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' rvFUN()

rvFUN <- function(x, obs, date, group, reflective = TRUE, date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00', theta = NA, tau = NA) {
  
  # install and load required packages
    list.of.packages <- c("stats", "stringr", "raster", "ggplot2", "data.table", "dplyr")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    library(stats); library(stringr); library(raster); library(ggplot2); library(data.table); library(dplyr)
  
    # clause on type of analysis to be run
    if (reflective == TRUE) {
      
      date.start = date.start 
      date.end = date.end
    } else {
      
      date.start = Sys.time()-60*60*24*7
      date.end = Sys.time()
    }
    
    date.interval <- interval(date.start, date.end, tzone = 'Pacific/Auckland')
    
    # define selected variables
    # use `data.table` package to deal with large datasets
    x <- as.data.table(x)
    x[, date := ymd_hms(date, tz = 'Pacific/Auckland')]
    x <- x[date %within% date.interval][, date := force_tz(date, tzone = 'Pacific/Auckland')]
    x$group <- x[,..group]
    x$obs <- x[,..obs]
	
	# add day column 
	x <- x[, day := str_sub(date, end = 10)]
  
	# add daily variance and daily number of group1 across day and group2
	x <- x[, day.var := var(obs, na.rm = T), by = list(group, day)]
	x <- x[, n.grp := length(unique(group)), by = list(day)]
  
	# select the columns for day, groups, daily variance & n of sites
	x <- x[, list(day, test = 'rank variance', group, day.var, n.grp)]
  
	# remove repeating and empty rows
	x <- unique.data.frame(x)
	x <- na.omit(x)
  
	# create ecdf column (this will create warnings)
	x <- suppressWarnings(x[, result := ecdf(day.var)(unique(day.var)), by = list(day)])
  
	# use theta and tau thresholds if present
	if((!is.na(theta))) {
  
	    # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
	    if (theta < 0.5) {
	      
			x <- x[, warning := ifelse(result < theta, 1, 0), by = list(group)]
	    
	    } else {
	      
			x <- x[, warning := ifelse(result > theta, 1, 0), by = list(group)]
	    
		}
	    
		if((!is.na(tau))) {
  
		# this looks at running means and is for tau, the length of time for alarms
	
			x <- x[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)] 
			 
		} else {
		
			x$alaram <- NA
		} 
		}  else {
	
			x$warning <- NA
			x$alarm <- NA
		}
  
	# gather variables of interest and return
	x <- x[, list(date = as.POSIXct(day), site = group, test = 'ranked variance', statistic = result, warning, alarm)]

	x

}
