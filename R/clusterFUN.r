#' A daily cluster function.
#'
#' This function gives the daily hierarchical cluster for a set of observations.
#' @param obs This is the column under observation.
#' @param date This the date column (set up as YYYY-MM-DD HH:MM:SS).
#' @param group This is the column where data will be partitioned upon (either site or time).
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' clusterFUN()

clusterFUN <- function(x, date, obs, group, reflective = TRUE, theta = NA, tau = NA, date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00'){
  
  list.of.packages <- c("plyr", "raster", "stringr", "lubridate", "data.table")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(plyr); library(raster); library(stringr); library(lubridate); library(data.table)
  
  # clause on type of analysis to be run, define dates
  if (reflective == TRUE) {
    
    date.start = date.start 
    date.end = date.end
  } else {
    
    date.start = Sys.time()-60*60*24*7
    date.end = Sys.time()
  }
  
  date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end)
  
  # define selected variables
  # use data.table package to deal with large datasets
  x <- as.data.table(x)
  setDT(x)
  x$date <- x[,..date]
  x[, date := ymd_hms(date)]
  x <- x[date >= date.start & date <= date.end]
  x$group <- x[,..group]; x$obs <- x[,..obs]
  
  # set up data
  cast.dat <- x[, list(date, obs, group)]
  cast.dat <- cast.dat[, dcast(cast.dat, ... ~ group, median, value.var = 'obs')]
  cast.dat <- setDT(cast.dat)[, day := str_sub(date, end = 10)]
 
  # run cluster function
  clusterA <- cast.dat %>%
	group_by(day) %>%
	filter(length(date) > 2) %>%
	do(distFUN(., date = 'date'))
	
  # use theta and tau thresholds if present
  if((!is.na(theta))) {
  
	clusterA <- setDT(clusterA)[, warning := ifelse(dissimilarity > theta, 1, 0), by = list(day)]
	  
	if((!is.na(tau))) {

		clusterA <- clusterA[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)]
				
	} else {
			
		clusterA$warning <- NA
	}
  } else {
		
	clusterA$warning <- NA; clusterA$alarm <- NA

  }
    
  # tidy up file and return
  clusterA <- clusterA[, list(date = as.POSIXct(day), test = 'cluster', site, statistic = dissimilarity, warning, alarm)]
  return(clusterA)

}