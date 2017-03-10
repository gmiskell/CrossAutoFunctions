#' A rolling Kolmogorov-Smirnov (KS) two-sample test function
#'
#' This function gives the rolling KS test for an observation, relative to a proxy from other observations.
#' @depen obs The time-series under investigation
#' @depen site The location of the time-series (or variable description)
#' @depen type This is used to identify the commonality in the network, and can allow the function to run over multiple networks
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param proxy This decides what type of proxy to use as the comparison. Options are 'network median' (default), or 'nearest site'.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @param window.length This defines the length of the sampled window, defined by row number.
#` @export
#' @examples
#' ksFUN()


ksFUN <- function(x, obs, site, type, reflective = TRUE, plot = TRUE, proxy = 'network median', date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00', theta = NA, tau = NA, window.length = 72){
  
  # install and load required packages
  list.of.packages <- c("raster", "ggplot2", "plyr", "data.table", "stringr", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  library(raster); library(ggplot2); library(plyr); library(data.table); library(stringr); library(dplyr)
    
  # clause on analysis type
  if (reflective == TRUE){
    plot.lab = 'reflective'
  } else {
    date.start = Sys.Date()-1
    date.end = Sys.Date()
    plot.lab = 'live'
  }
  
  # select required data
  a <- x %>%
	  select(date, site, obs, type) %>%
	  filter(date >= ymd_hms(date.start), date <= ymd_hms(date.end)) %>%
    arrange(date) %>%
    na.omit()
	
  # run selected proxy function (either network median or nearest site)
  if(proxy == 'network median'){
    b <- as.data.table(a)
    b <- b[, proxy := networkmedianFUN(b), by = list(date, type)]
	  proxy_selection = 'network median'
  }
  if(proxy == 'nearest site'){
	  b <- as.data.table(a)
    b <- b[, proxy := nearestsiteFUN(b), by = list(date, type)]
    proxy_selection = 'nearest site'
  }
  
  # generate rolling KS test results using function by site
  c <- ddply(b, .(site),
                    function(x) rollingKStest(x, window = window.length))
  
  # filter to those with p-values < 0.05 (significantly different to the network proxy)
  c$theta.alarm <- ifelse(c$p.value < theta, 1, 0)
  
  #alarm generated is for systematic or local change between obs and proxy using window.2
  d <- ddply(c, .(site),
                  transform,
                  tau.alarm = movingFun(theta.alarm, n = tau, fun = 'mean', type = 'to'))
  
  # gather variables of interest and return
  ks.filter <- d %>%
    mutate(day = str_sub(date, end = -9L), test = 'ks filter') %>%
    select(date, obs, site, p.value, test, filter) %>%
    unique()
  
  ks.filter <- cbind(ks.filter, proxy_selection, theta, tau)
  
  if(plot == TRUE) {
    plot <- ggplot(ks.filter, aes(x = date, y = p.value)) +
      geom_line(aes(col = type)) +
	    facet_wrap(~site)
	  return(c(plot, ks.filter))
  }
  
}
  