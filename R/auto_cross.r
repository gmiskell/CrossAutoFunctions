#' A function to set up data ready to perform "auto" analysis
#'
#' This function reformats the data relative to auto analysis, where site is held steady.
#' @param x This is the column under observation.
#' @param group1 This is the id of each set of data being compared.
#' @param ell This sets the length of time backwards (in days) that should be considered as the auto proxy. Defaults to 3.
#' @param group2 This is the column which identifies separate networks within the data.
#' @param date.start This sets the start of the selected data. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data. Defaults to one week later than `date.start`.
#' @export
#' @examples
#' autoFUN()

autoFUN <- function(x, ell = 3, group1 = 'site', group2 = 'pol', date.start = '2016-07-01', date.end = NA){
  
  list.of.packages <- c("plyr", "raster", "stringr", "data.table", "lubridate", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
 
  library(plyr); library(raster); library(data.table); library(stringr); library(lubridate); library(dplyr)
  
  date.start.ell <- as.character(as.Date(date.start) - days(ell))
  if(is.na(date.end) == TRUE) {date.end <- as.character(as.Date(date.start) + days(7))}
  
  x$group1 <- x[, ...group1]; x$group2 <- x[, ...group2]
  
  x <- x %>%
	  select(date, group1, group2, value) %>%
    filter(date >= (date.start.ell), date < (date.end)) %>%
    group_by(group1) %>%
	  filter(!is.na(group2)) %>%
	  arrange(date)
	len <- as.numeric(difftime(x$date[2], x$date[1], units = 'mins'))
	if(len == 60) len = 24
	if(len == 1) len = 1440	
	nn <- seq(1:group1) * len		
	x <- setDT(x)[, paste("l", 1:ell) := shift(value, n =  nn, type = 'lag'), by = group1]
	x <- data.frame(x)
	
  y <- x %>%
	  rename(response.value = value) %>%
	  group_by(date, group1) %>%
	  melt(id.vars = c('date', 'group1', 'group2', 'response.value')) %>%
	  rename(comparison = variable, comparison.value = value)
  
  y

}

#' A function to set up data ready to perform "cross" analysis
#'
#' This function reformats the data relative to cross analysis, where time is held steady. NB: id & group1 here are the same, and so only id is included.
#' @param x This is the column under observation.
#' @param group1 This is the id of each set of data being compared.
#' @param group2 This is the column which identifies separate networks within the data.
#' @param date.start This sets the start of the selected data. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data. Defaults to one week later than `date.start`.
#' @export
#' @examples
#' crossFUN()

crossFUN <- function(x, group1, group2, date.start = '2016-07-01', date.end = NA){
  
  list.of.packages <- c("reshape2", "lubridate", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
 
  library(reshape2); library(lubridate); library(dplyr)
  
  if(is.na(date.end) == TRUE) {date.end <- as.character(as.Date(date.start) + days(7))}
  
  x$group1 <- x[, ...group1]; x$group2 <- x[, ...group2]
  
  x <- x %>%
    select(date, group1, value, group2) %>%
  	filter(date >= (date.start) & date <= (date.end))

  cast.x <- dcast(x, ... ~ group1)
  
  y <- x %>%
    join(cast.x, by = c('date', 'group2')) %>%
    rename(response.value = value) %>%
    group_by(date, group1) %>%
    melt(id.vars = c('date', 'group1', 'group2', 'response.value')) %>%
    rename(comparison = variable, comparison.value = value)
  
  y <- y[!(y$group1 == y$comparison),]
	
	y
	
}
