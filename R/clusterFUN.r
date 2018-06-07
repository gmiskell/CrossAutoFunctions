#' A daily cluster function.
#'
#' @param obs The column under observation.
#' @param date A date column (set up as `YYYY-MM-DD HH:MM:SS`).
#' @param group A column which the data are partitioned on (either spatial or temporal).
#' @param reflective Decides whether to use a selected time period of data (reflective, `TRUE`, default), or to use data from the latest day (live, `FALSE`).
#' @param by.day An option where analysis uses data from only that day ("cross", `TRUE`, default), or uses all selected data ("auto", `FALSE`). If selecting auto, the last five days will be used.
#' @param theta Sets the test threshold on whether to flag the data. Defaults to `NA` (no flags given).
#' @param tau Sets the day threshold on whether to flag the data given consistent theta flags. Defaults to `NA` (no flags given).
#' @export
#' @examples
#' clusterFUN()

clusterFUN <- function(x, date, obs, group, by.day = TRUE, reflective = TRUE, theta = NA, tau = NA){
  
  library(raster);library(tidyverse);library(lubridate);library(data.table);
     
  # define selected variables
  # use data.table package to deal with large datasets
  x <- as.data.table(x);
  setDT(x);
  x$date <- x[[date]];
  x[, date := ymd_hms(date)];
  
  # clause on type of analysis to be run, define dates
  if(reflective == TRUE) {
    date.start = min(x$date, na.rm = T); 
    date.end = max(x$date, na.rm = T);
  } else {  
    date.start = now()-60*60*24*7;
    date.end = now();
    date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end);
  };
  
  x <- x[date %within% interval(date.start, date.end)];
  x$group <- x[[group]]; x$obs <- x[[obs]];
  
  if(by.day == TRUE){
    # set up data
    cast.dat <- x[, list(date, obs, group)];
    cast.dat <- unique(cast.dat);
    cast.dat <- spread(cast.dat, obs, group, -date);
    cast.dat$day <- str_sub(cast.dat$date, end = 10);
 
    # run cluster function
    clusterA <- cast.dat %>%
	  group_by(day) %>%
	  filter(length(date) > 2) %>%
	  do(distFUN(.));
  
    # use theta and tau thresholds if present
    if(!is.na(theta)) {
      clusterA <- setDT(clusterA)[, warning := ifelse(dissimilarity > theta, 1, 0)];
	  if(!is.na(tau)) {
        clusterA <- clusterA[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)]; 
      } else {  
        clusterA$warning <- NA
      };
    } else {
      clusterA$warning <- NA; clusterA$alarm <- NA;
    };
  
    # tidy up file and return
    clusterA <- clusterA[, list(date = as.POSIXct(day), test = 'cluster', group, statistic = dissimilarity, warning, alarm)];
  };
  
  if(by.day == FALSE){
    # set up data
    cast.dat <- x[, list(date, obs, group)];
    cast.dat <- unique();
    cast.dat <- spread(cast.dat, obs, group, -date);
    cast.dat$day <- str_sub(cast.dat$date, end = 10);
    
    # run cluster function
    clusterA <- cast.dat %>%
      dplyr::select(-date) %>%
      nest(-day) %>%
      mutate(num.col = map_int(data, ~sum(map_lgl(.x, ~!all(is.na(.x)))))) %>%
      unnest() %>%
      filter(num.col > 2 & length(day) > 2);
    
    if(length(clusterA$day) > 0){
      clusterA <- clusterA %>%
      dplyr::select(-num.col) %>%
      group_by(day) %>%
      do(distFUN(.)) %>%
      filter(group == 'l 0');
    
    # use theta and tau thresholds if present
    if(!is.na(theta)){
	  clusterA <- setDT(clusterA)[, warning := ifelse(dissimilarity > theta, 1, 0)];
	  if(!is.na(tau)) {
	 	clusterA <- clusterA[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)];		
	  } else {	
		clusterA$warning <- NA;
	  }
  } else {
    clusterA$warning <- NA; clusterA$alarm <- NA;
  }
    
  # tidy up file and return
  clusterA <- clusterA[, list(date = as.POSIXct(day), test = 'cluster', group, statistic = dissimilarity, warning, alarm)];
  } else {
    clusterA <- data.frame(date=as.Date(character()),
                                 test=character(),
                                 group=character(),
                                 statistic=numeric(),
                                 warning=numeric(),
                                 alarm=numeric(),
                                 stringsAsFactors=FALSE);
  }
}
    return(clusterA);
}
