#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#` @export
#' @examples
#' rvFUN()

rvFUN <- function(x, obs, reflective = TRUE, plot = TRUE, date.start = '2016-07-01', date.end = '2016-07-10', theta = NA, tau = NA){
  
  # load required libraries
  library(stats); library(stringr); library(raster); library(data.table); library(dplyr)
  
  # clause on type of analysis to be run
  if (reflective == TRUE) {
  date.start = date.start 
  date.end = date.end
  plot.lab = 'reflective'
  } else {
  date.start = str_c(Sys.Date()-1)
  date.end = str_c(Sys.Date())
  plot.lab = 'live'
  }
  
  # this is a package in R that helps deal with large datasets
  x <- as.data.table(x)
  
  # add day column 
  setDT(x)[, day := str_sub(date, end = 10)]
  
  # add daily variance, daily number of sites  and ecdf column (by site and pol type)
  x <- x[, day.var := var(-conc_mean, na.rm=T), by = list(site, pol, day)]
  x <- x[, n.site := length(unique(site)), by = list(day, pol)]
  
  # select the columns for day, site, PM type, daily variance & n of sites
  x <- x[, list(day, test, site, pol, day.var, n.site)]
  
  # remove repeating and empty rows
  x <- unique.data.frame(x)
  x <- na.omit(x)
  
  # create ecdf column (this will create warnings)
  x <- suppressWarnings(x[, result := ecdf(day.var)(unique(day.var)), by = list(day, pol)])
  
  # use theta and tau thresholds
  if((!is.na(theta))) {
  
  if((!is.na(tau))) {
  
    if (theta < 0.5){
    
    # return back to df to use a function
    x <- x %>%
    as.data.frame() %>%
    # put data into distinct site and PM type groups
    group_by(site, pol) %>%
    # if else clause on whether the data are less than theta, becomes 1 if true, and 0 otherwise
    mutate(change.thresh = ifelse(result < theta, 1, 0)) %>%
    # this looks at running means and is for tau, the length of time for alarms
    mutate(alarm = movingFun(change.thresh, n = tau, type = 'to', fun = mean, na.rm = T)) 
  } else {
    
    # this is similar to when theta < 0.5
    x <- x %>%
    as.data.frame() %>%
    group_by(site, pol) %>%
    mutate(change.thresh = ifelse(result > theta, 1, 0)) %>%
    mutate(alarm = movingFun(change.thresh, n = tau, type = 'to', fun = mean, na.rm = T)) 
  }
  } else {
  tau = NA
  }
  }  else {
  theta = NA
  tau = NA
  }
  
  ggplot(x, aes(x = as.POSIXct(day), y = result, colour = site)) +
	geom_line() +
	facet_wrap(~site) +
	
  # returns z as product of function
  x

}
    

    
  
      
      