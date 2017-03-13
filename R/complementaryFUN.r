#' A test to check if values are too high.
#'
#' This function gives binary output on whether values are over a defined threshold.
#' @param value This is the column of interest, with default being 'value'.
#' @param theta This is the threshold of what a 'high' value is.
#` @export
#' @examples
#' hcFUN()

# complementary functions for Christchurch's network analysis

# high concentration script
hcFUN <- function(x, value = value, theta = 50){
  
  library(data.table)
  
  x <- data.table(x)
  x[, test := 'high concentration', filter := ifelse(value > theta, 1, 0)]
  
}

#' A test to check where the flow rate is too low
#'
#' @param theta This is the threshold of what 'low' flow is.
#` @export
#' @examples
#' frFUN()

# low flow rate script

frFUN <- function(x, theta = 1.5){
  
  library(data.table)
  x <- data.table(x)
  x <- x[, filter := ifelse(flow < theta | is.na(flow), 1, 0), test := 'flow']
 
}

#' A test to see how long it has been since visits.
#'
#' @param theta This is the threshold of when length between site visits is 'high' (in days).
#` @export
#' @examples
#' vlFUN()


vlFUN <- function(x, theta = 60){

	library(lubridate); library(dplyr)
	
	x <- x %>%
		group_by(site, pol) %>%
		mutate(length = difftime(date, last.visit, units = 'hours'), test = 'last visit') %>%
		mutate(filter = ifelse(length > length.threshold, 1, 0))
		
}  

