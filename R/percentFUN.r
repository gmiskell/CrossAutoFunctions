#' A function for comparing the daily minima and maxima.
#'
#' This function finds and compares the daily 'minima' (5th) and 'maxima' (95th) percentiles across groups using modified z-scores (Mean Absolute Deviation - MAD).
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param group This is the first grouping value of the data, with the default being 'site'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' #' @param by.day This is an option for where analysis uses data from only that day ('cross', default), or uses all selected data ('auto', FALSE). If selecting auto, the last five days will be used.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to 2 (similar to 2 standard deviations).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' export
#' examples
#' percentFUN()

percentFUN <- function(x, date, obs, group, by.day = T, reflective = TRUE, theta = 2, tau = NA){
  
 library(data.table);library(tidyverse);
  
  # define selected variables
  # use `data.table` package to deal with large datasets
  x <- as.data.table(x);
  x$date <- x[[date]];
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
	
  # determine the 5th and 95th percentiles
  min.test <- x[, quart := quantile(obs, probs = 0.05, na.rm = T), by = list(day, group)][, list(day, test = 'min_percentile', group, quart)];
  max.test <- x[, quart := quantile(obs, probs = 0.95, na.rm = T), by = list(day, group)][, list(day, test =  'max_percentile', group, quart)];
  min.max <- rbind(min.test, max.test);

  # remove repeating rows
  min.max <- unique(min.max);

  # derive daily MAD scores
  min.max <- min.max[, median := median(quart, na.rm = T), by = list(day, test)][, sd := mad(quart, na.rm = T), by = list(day, test)];
  
  # determine MAD score effect sizes
  min.max.test <- min.max[, sigma := (quart - median)/(1.4826 * (sd+0.0000001))];
  };
  
  if(by.day == FALSE){
    # add day column 
    x <- x[, day := str_sub(date, end = 10)];
    
    # determine the 5th and 95th percentiles
    min.test <- x[, quart := quantile(obs, probs = 0.05, na.rm = T), by = list(day)][, list(day, test = 'min_percentile', quart)];
    max.test <- x[, quart := quantile(obs, probs = 0.95, na.rm = T), by = list(day)][, list(day, test =  'max_percentile', quart)];
    min.max <- rbind(min.test, max.test);
    
    # remove repeating rows
    min.max <- unique(min.max);
    
    # derive daily MAD scores - running over five days
    min.max <- min.max[, med := movingFun(quart, 5, median, na.rm = T, type = 'to')][, sd := movingFun(quart, 5, mad, na.rm = T, type = 'to')];
    
    # determine MAD score effect sizes
    min.max.test <- min.max[, sigma := (quart - med)/(1.4826 * (sd+0.0000001))];
  }
  # use theta and tau thresholds
  if(!is.na(theta)){
  
    min.max.test <- min.max.test[, warning := ifelse(abs(sigma) > theta, 1, 0)];
    
	if(!is.na(tau)){
  
    # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
    # this looks at running means and is for tau, the length of time for alarms
		min.max.test <- min.max.test[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T), by = list(test)];
	} else {
	
		min.max.test$alarm <- NA;
	}; 
  } else {
	
    min.max.test$warning <- NA;
	  min.max.test$alarm <- NA;
  };
  
  # gather variables of interest and return
  min.max.test <- min.max.test[, list(date = as.POSIXct(day), group, test, statistic = sigma, warning, alarm)];
  return(min.max.test);

};
