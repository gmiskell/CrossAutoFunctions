#' A daily cluster function
#'
#' This function gives the daily hierarchical cluster for a set of observations. Require date column (as.POSIXct), site column (character), pol column (character), value column (numeric)
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#` @export
#' @examples
#' clusterFUN()


clusterFUN <- function(x, obs, reflective = TRUE, plot = TRUE, theta = NA, tau = NA, date.start = '2016-07-01', date.end = '2016-07-10'){
  
  list.of.packages <- c("plyr", "raster", "stringr", "lubridate", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
 
  library(plyr); library(raster); library(stringr); library(lubridate); library(dplyr)
  
  # clause on type of analysis to be run
  if (reflective == TRUE){
    date.start = date.start
    date.end = date.end
    date.start = with_tz(as.POSIXct(date.start), tzone = time.zone)
    date.end = with_tz(as.POSIXct(date.end), tzone = time.zone)
	plot.lab = 'reflective'
  } else {
    date.start = str_c(Sys.Date()-1)
    date.end = str_c(Sys.Date())
	date.start = with_tz(as.POSIXct(date.start), tzone = time.zone)
    date.end = with_tz(as.POSIXct(date.end), tzone = time.zone)	
	plot.lab = 'live'
  }
  
  # set up data
  cast.dat <- x %>%
    select(date, obs, site, pol) %>%
    dcast(... ~ site + pol, median)
  cast.dat$day <- day(cast.dat$date)	
  
  # generate output using distFUN
  if(plot == TRUE){ 
    # turn on pdf generator so to combine all plots
    plot.name <- str_c('cluster_plot_', plot.lab, '.pdf')
    pdf(plot.name, width = 6, height = 10)
  
    # apply the function within the function
    if(auto == TRUE) {
	d_ply(cast.dat, .(day, site),
					  function(x) distFUN(x, plot = TRUE))

	} else {
	d_ply(cast.dat, .(day, pol),
                      function(x) distFUN(x, plot = TRUE))  
    }					  
    
	dev.off()
	dev.off()
	
  } else {
    if(auto == TRUE) {
		clusterAuto <- cast.dat %>%
		group_by(day, site) %>%
		do(distFUN(., plot = FALSE))
	if(!is.na(tau)) {
	clusterAuto <- clusterAuto %>%
	    group_by(day, group) %>%
		select(site, group) %>%
		unique() %>%
		mutate(n = n()) %>%
		ungroup() %>%
		mutate(alarm = ifelse(n <= 1, 1, 0)) %>%
		group_by(site) %>%
		mutate(filter = movingFun(alarm, n = tau, type = 'to', fun = mean)) 
	}		
    write.csv(clusterAuto, filename = str_c('AUTO_cluster', plot.lab, '.csv'))
	} else {	
	    clusterCross <- cast.dat %>%
	    group_by(day, pol) %>%
	    do(distFUN(., plot = FALSE))	  
    write.csv(clusterCross, filename = str_c('CROSS_cluster_', plot.lab, '.csv'))
  }
}   
}
  
  
  
  