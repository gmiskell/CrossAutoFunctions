#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' hcFUN()

# complementary functions for Christchurch's network analysis

# high concentration script
hcFUN <- function(x, value = conc_mean, high.conc.threshold = 50){
  
  library(data.table)
  
  x <- data.table(x)
  x[, test := 'high concentration', filter := ifelse(value > high.conc.threshold, 1, 0)]
  
  x <- x %>%
    mutate(test = 'high concentration', filter = ifelse(value > high.conc.threshold, 1, 0))
     
}

#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' frFUN()

# low flow rate script

frFUN <- function(x, flow.threshold = 1.5){
  
  library(data.table)
  
  x <- x %>%
    mutate(filter = ifelse(flow < flow.threshold | is.na(flow), 1, 0), test = 'flow') 
 
}

#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' vlFUN()

# visit length script
# length of time until 'too long' between site visits (in days)

vlFUN <- function(x, length.threshold = 60){

	library(lubridate); library(dplyr)
	
	x <- x %>%
		group_by(site, pol) %>%
		mutate(length = difftime(date, last.visit, units = 'hours'), test = 'last visit') %>%
		mutate(filter = ifelse(length > length.threshold, 1, 0))
		
}

#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' metFUN()
	
# met filter

metFUN <- function(x, temp.threshold = 2, ws.threshold = 1){
  
  library(dplyr)
   
  x <- x %>%
    mutate(filter = ifelse(temp < temp.thresh | ws < ws.thresh, 1, 0),
           test = 'meteorology')
  
}
  
  

