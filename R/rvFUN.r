#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param group This is the first grouping value of the data, with the default being 'site'.
#' @param reflective This decides whether to use all selected data (TRUE), or to use data from the latest day (live). Defaults to TRUE.
#' @param by.day This is an option for where analysis uses data from only that day ('cross', default), or uses all selected data ('auto', FALSE). If selecting auto, the last five days will be used.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' rvFUN()

rvFUN <- function(x, obs, date, group, ell = 5, reflective = TRUE, by.day = TRUE, theta = NA, tau = NA) {
  
library(zoo);library(data.table);library(tidyverse);
    
    # define selected variables
    # use `data.table` package to deal with large datasets
    x <- as.data.table(x);
    x[, date := ymd_hms(date)];
    
     # clause on type of analysis to be run
    if (reflective == TRUE){
      
      date.start = min(x$date, na.rm = T);
      date.end = max(x$date, na.rm = T);
    } else {
      
      date.start = Sys.time()-60*60*24*7;
      date.end = Sys.time();
      date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end);
    
    };
    
    x <- x[date %within% interval(date.start, date.end)];
    x$group <- x[[group]];
    x$obs <- x[[obs]];
	
  if(by.day == TRUE){
    
	  # add day column
	  x <- x[, day := str_sub(date, end = 10)];
  
	# add daily variance and daily number of group1 across day and group2
	x <- x[, day.var := var(obs, na.rm = T), by = list(group, day)];
	x <- x[, n.grp := length(unique(group)), by = list(day)];
  
	# select the columns for day, groups, daily variance & n of sites
	x <- x[, list(day, test = 'rank variance', group, day.var, n.grp)];
  
	# remove repeating and empty rows
	x <- unique.data.frame(x);
	x <- na.omit(x);
  
	# create ecdf column (this will create warnings)
	x <- suppressWarnings(x[, result := ecdf(day.var)(unique(day.var)), by = list(day)]);
	
  };
  
  if(by.day == FALSE){
    
    x <- x[, day := str_sub(date, end = 10)];
    
    # add grouped variance and daily number of group1 across day and group2
    x <- x[, day.var := var(obs, na.rm = T), by = list(group)];
    
    # select the columns for day, groups, daily variance & n of sites
    x <- x[, list(day, test = 'rank variance', group, day.var)];
    
    # remove repeating and empty rows
    x <- unique.data.frame(x);
    x <- na.omit(x);
    
    # create ecdf column (this will create warnings)
    x.zoo <- zoo(x$day.var);
    ecdf <- rollapply(x.zoo, function(x.zoo) ecdf(x.zoo)(unique(x.zoo)), 
                      by.column = F, fill = NA, align = 'right', width = ell);
    x$result <- as.numeric(ecdf[, ncol(ecdf)]);
    
  }
    
  
	# use theta and tau thresholds if present
	if((!is.na(theta))){
  
	    # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
	    if (theta < 0.5){
	      
			x <- x[, warning := ifelse(result < theta, 1, 0)];
	    
	    } else {
	      
			x <- x[, warning := ifelse(result > theta, 1, 0)];
	    
		};
	    
		if((!is.na(tau))){
  
		# this looks at running means and is for tau, the length of time for alarms
	
			x <- x[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)]; 
			 
		} else {
		
			x$alaram <- NA;
		};
		}  else {
	
			x$warning <- NA;
			x$alarm <- NA;
		};
  
	# gather variables of interest and return
	x <- x[, list(date = as.POSIXct(day), group, test = 'ranked variance', statistic = result, warning, alarm)];

	x;

};
