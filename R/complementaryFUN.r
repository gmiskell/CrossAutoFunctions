#' A test to check if values are too high.
#'
#' This function gives binary output on whether values are over a defined threshold.
#' @param obs This is the column of interest, with default being 'value'.
#' @param theta This is the threshold of what a 'high' value is.
#' @export
#' @examples
#' hcFUN()

hcFUN <- function(x, obs = 'value', theta = 50){
  
  list.of.packages <- c("data.table")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(data.table)
  
  x <- as.data.table(x)
  x$obs <- x[,..obs]
  x <- x[, test := 'high concentration'][, alarm := ifelse(obs > theta, 1, 0)]
  
  x <- x[, list(date, test, alarm)]
  
}

#' A test to check where the flow rate is too low
#'
#' @param theta This is the threshold of what 'low' flow is.
#' @param flow This is the flow column.
#' @export
#' @examples
#' frFUN()

frFUN <- function(x, flow = 'flow', theta = 1.5){
  
  list.of.packages <- c("data.table")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(data.table)
  
  x <- as.data.table(x)
  x$flow <- x[,..flow]
  x <- x[, alarm := ifelse(flow < theta | is.na(flow), 1, 0)][, test := 'flow']
 
  x <- x[, list(date, test, alarm)]
  
}

#' A test to see how long it has been since visits.
#'
#' @param theta This is the threshold of when length between site visits is 'high' (in days).
#' @param date This is the date column (in YYYY-MM-DD).
#' @param last.visit This is the column of the date of last visit to the site (in YYYY-MM-DD).
#' @export
#' @examples
#' vlFUN()

vlFUN <- function(x, date, last.visit = 'last.visit', theta = 60){

  list.of.packages <- c("lubridate", "data.table")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(lubridate); library(dplyr)
	
  x <- as.data.table(x)
  x$date <- x[,..date]; x$last.visit <- x[,..last.visit]
  x <- x[, test := 'last visit'][, alarm := ifelse(difftime(date, last.visit, unit = 'days') > theta, 1, 0)]
  
  x <- x[, list(date, test, alarm)]
		
}  
