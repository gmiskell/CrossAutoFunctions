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
#' clusterFUN()

####		rolling cluster function script

# description:	this filter will compare daily clusters of locations
#              	using hierarchical cluster methods. Options to make this run on live data
#				or on past data (reflective) and to make tests auto (use data from one
#				location) or cross (use data from multiple locations)

# data set-up:	require - date column (as.POSIXct), site column (character), pol
#				column (character), value column (numeric)

# options:		x - data
#				reflective =	TRUE (default), use all data within date.start and
#								date.end span, or FALSE use data from latest day only
#				plot = 			TRUE (default), will create plots, or FALSE will create
#								csv file
#				theta = 		0.3 (1-theta) to give dissimilarity threshold
#				tau =			3 (days) to give consistency threshold
#				date.start =	'2016-01-01 00:00:00'
#				date.end =		'2017-01-01 00:00:00'

clusterFUN <- function(x, reflective = TRUE, plot = TRUE, theta = 0.3,
						tau = 3, ell = 5, date.start = '2016-01-01 00:00:00', date.end = '2017-01-01 00:00:00', time.zone = Sys.timezone()){
  
  # load required libraries
  library(raster); library(stringr); library(lubridate); library(dplyr)
  
  # clause on type of analysis to be run
  if (reflective == TRUE){
    date.start = date.start
    date.end = date.end
    date.start = with_tz(as.POSIXct(date.start), tzone = time.zone)
    date.end = with_tz(as.POSIXct(date.end), tzone = time.zone)
	plot.lab = 'reflective'
  } else {
    date.start = str_c(Sys.Date()-1, ' 00:00:00')
    date.end = str_c(Sys.Date(), ' 00:00:00')
	date.start = with_tz(as.POSIXct(date.start), tzone = time.zone)
    date.end = with_tz(as.POSIXct(date.end), tzone = time.zone)	
	plot.lab = 'live'
  }
  
  # clause on either auto or cross test to be run
  if (auto == TRUE){
  cast.dat <- x %>%
    filter(date >= as.POSIXct(date.start) - days(ell), date < (date.end)) %>%
	mutate(day = str_sub(date, end = 10)) %>%
    group_by(site) %>%
	filter(!is.na(pol)) 	
	len <- as.numeric(difftime(cast.dat$date[2],cast.dat$date[1],units = 'mins'))
	if(len == 60) len = 24
	if(len == 1) len = 1440	
	nn <- seq(1:ell) * len		
	cast.dat <- setDT(all.data.melt.hr)[, paste("l", 1:ell) := shift(value, n =  nn),by = site]
	cast.dat <- data.frame(cast.dat)
  } else {
  # rearrange data structure and create day variable
  cast.dat <- x %>%  
    filter(date >= (date.start) & date < (date.end)) %>%
    dcast(... ~ site, value.var = 'value', median) %>%
    mutate(day = str_sub(date, end = 10L))
  }

  # generate output using distFUN function
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
  
  
  
  