#' A daily cluster function.
#'
#' @param obs This is the column under observation.
#' @param date This the date column (set up as YYYY-MM-DD HH:MM:SS).
#' @param group This is the column where data will be partitioned upon (either site or time).
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param by.day This is an option for where analysis uses data from only that day ('cross', default), or uses all selected data ('auto', FALSE). If selecting auto, the last five days will be used.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' clusterFUN()

clusterFUN <- function(x, date, obs, group, by.day = TRUE, reflective = TRUE, theta = NA, tau = NA){
  
  list.of.packages <- c("plyr","purrr","liraster","stringr","lubridate","data.table");
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])];
  if(length(new.packages)) install.packages(new.packages);
  library(plyr);library(purrr);library(raster);library(stringr);library(lubridate);library(data.table);
     
  # define selected variables
  # use data.table package to deal with large datasets
  x <- as.data.table(x);
  setDT(x);
  x$date <- x[,..date];
  x[, date := ymd_hms(date)];
  
  # clause on type of analysis to be run, define dates
  if (reflective == TRUE) {
    
    date.start = min(x$date, na.rm = T); 
    date.end = max(x$date, na.rm = T);
  } else {
    
    date.start = Sys.time()-60*60*24*7;
    date.end = Sys.time();
    date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end);
  };
  
  x <- x[date %within% interval(date.start, date.end)];
  x$group <- x[,..group]; x$obs <- x[,..obs];
  
  if(by.day == TRUE){
    # set up data
    cast.dat <- x[, list(date, obs, group)];
    cast.dat <- cast.dat[, dcast(cast.dat, ... ~ group, median, value.var = 'obs')];
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
    cast.dat <- cast.dat[, dcast(cast.dat, ... ~ group, median, na.rm = T, value.var = 'obs')];
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
