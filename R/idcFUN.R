#' An inverse distance correlation function.
#' 
#' This function compares the correlations of one site against the other sites in the network. The final product is a sum of the correlations multiplied by the inverse distance. Data are rolled iteratively to make a continuous examination of network correlations over time.
#' @param obs The column under examination.
#' @param date A date column (set up as YYYY-MM-DD HH:MM:SS).
#' @param group This is the column where data will be partitioned upon (either site or time).
#' @param lat Latitude.
#' @param lon Longitdue.
#' @param reflective This decides whether to use a selected time period of data (reflective, `TRUE`, the default setting), or to use the data from the latest day (live, `FALSE`).
#' @param theta This sets the test threshold on whether to flag the data. Defaults to `NA` (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to `NA` (no flags given).
#' @export
#' @examples 
#' idcFUN()


idcFUN <- function(x, obs = 'obs', date = 'date', group, lat, lon, reflective = TRUE, theta = 0.5, tau = 120, window.size = 72){
  
  list.of.packages <- c("raster","zoo","tidyverse","lubridate","data.table");
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])];
  if(length(new.packages)) install.packages(new.packages);
  library(raster);library(zoo);library(tidyverse);library(lubridate);library(data.table);
  
  # define selected variables
  # use data.table package to deal with large datasets
  x <- as.data.table(x);
  setDT(x);
  x$date <- x[,..date];
  x[, date := ymd_hms(date)];
  
  # clause on type of analysis to be run, define dates
  if(reflective == TRUE){
    date.start = min(x$date, na.rm = T);
    date.end = max(x$date, na.rm = T);
    } else {
      date.start = Sys.time()-60*60*24*7;
      date.end = Sys.time();
      date.start <- ymd_hms(date.start); date.end <- ymd_hms(date.end);
      };
  
  x <- x[date %within% interval(date.start, date.end)];
  x$group <- x[,..group]; x$obs <- x[,..obs]; x$lat <- x[,..lat]; x$lon <- x[..lon];
  
  spread.x <- x %>%
    dplyr::select(date, group, obs) %>%
    spread(obs, group, -date);
  
  date.x <- spread.x$date;
  date.x <- date.x[-c(1:(theta-1))];
  
  z <- zoo(spread.x[, -date]);
  c <- cor(z, use = 'pairwise.complete.obs');
  ut <- as.logical(c);
  n <- paste(rownames(c)[row(c)[ut]], rownames(c)[col(c)[ut]]);
  
  roll.corFUN <- rollapply(z, width = window.size, function(Z) {return(cor(Z)[ut])}, by.column = F, align = 'right');
  colnames(roll.corFUN) <- n;
  roll.corFUN.df <- data.frame(date.x, roll.corFUN);
  
  distances <- as.matrix(dist(cbind(x$lon, x$lat)));
  distance <- data.frame(distances);
  distance$group <- x$group;
  distance <- distance[!(is.na(distance$group) | distance$group==""),];
  distance <- Filter(function(x)!all(is.na(x)), distance);
  colnames(distance)[1:length(x$group)] <- as.character(distance$group);
  dist.melt <- distance %>%
    gather(group, dist) %>%
    mutate(distance = dist * 100);
  
  roll.corFUN.melt <- roll.corFUN.df %>%
    gather(date.x, correlation) %>%
    filter(correlation < 1) %>%
    left_join(dist.melt) %>%
    filter(distance > 1) %>%
    mutate(filt = correlation * 1/distance);
  
  # use theta and tau thresholds if present
  if(!is.na(theta)){
    roll.corFUN.melt <- setDT(roll.corFUN.melt)[, warning := ifelse(filt < theta, 1, 0)];
    if(!is.na(tau)){
      roll.corFUN.melt <- roll.corFUN.melt[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)];	
      } else {	
        roll.corFUN.melt$warning <- NA;
        };
    } else {
      roll.corFUN.melt$warning <- NA; roll.corFUN.melt$alarm <- NA;
      };
  
  # tidy up file and return
  roll.corFUN.melt <- roll.corFUN.melt[, list(date = as.POSIXct(date.x), test = 'idc', group, statistic = filt, warning, alarm)];
}
