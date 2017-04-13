#' A daily cluster function
#'
#' This function gives the daily hierarchical cluster for a set of observations. Require date column (as.POSIXct), grouping column (character - to be defined in the `autoFUN` function), second grouping column (character - if there are two or more networks operating at once), id column (character), and a value column (numeric).
#' @param x This is the column under observation.
#' @param group1 This is the column where data will be partitioned upon (either site or time).
#' @param group2 This is the column which identifies separate networks within the data.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' clusterFUN()


clusterFUN <- function(x, obs, id, group1, group2, reflective = TRUE, plot = TRUE, theta = NA, tau = NA, date.start = '2016-07-01', date.end = '2016-07-10'){
  
  list.of.packages <- c("plyr", "raster", "stringr", "lubridate", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(plyr); library(raster); library(stringr); library(lubridate); library(dplyr)
  
  # clause on type of analysis to be run
	if (reflective == TRUE) {
	
		date.start = date.start 
		date.end = date.end
	} else {
	
		date.start = str_c(Sys.Date()-7)
		date.end = str_c(Sys.Date())
	}
  
	# set up data
		cast.dat <- x %>%
			select(date, obs, group1, group2) %>%
			dcast(... ~ group1 + group2, median)
		cast.dat$day <- day(cast.dat$date)	
  
	
		clusterA <- cast.dat %>%
			group_by(day, group1) %>%
			do(distFUN(., d.thresh = theta, plot = FALSE))
	
	# use theta and tau thresholds if present
		if((!is.na(theta))) {
  
			if((!is.na(tau))) {

				clusterA <- clusterA %>%
					group_by(day, group1) %>%
					unique() %>%
					mutate(n = n()) %>%
					ungroup() %>%
					mutate(alarm = ifelse(n <= theta, 1, 0)) %>%
					group_by(group1) %>%
					mutate(filter = movingFun(alarm, n = tau, type = 'to', fun = mean)) 
			} else {
			
			tau <- NA
			}
		} else {
		
		theta <- NA
		tau <- NA
		}
    
	csv <- write.csv(clusterA, file = str_c(dest, '/cluster', start.date, '.csv'))
  
	# create plot and save
  
	# generate output using distFUN
		if(plot == TRUE){ 
    
		# turn on pdf generator so to combine all plots
			plot.name <- str_c('/cluster', start.date, '.pdf')
			pdf(plot.name, width = 6, height = 10)
  
		# apply the function within the function
			ts <- d_ply(cast.dat, .(day, id),
							function(x) distFUN(x, d.thresh = theta, plot = TRUE))
				      
			dev.off()
			return(c(ts, csv))
		} else { 
  
		# returns x as saved csv file
			csv
  }
}