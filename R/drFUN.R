#' A decomposed regression analysis
#'
#' This function gives the decomposed time series with long-term (week) and short-term (day) trends removed.
#' @param obs This is the numeric column under investigation.
#' @param date This is a required column with the format 'YYYY-MM-DD'.
#' @param group1 This is the first grouping value of the data, with the default being 'site'.
#' @param group2 This is the second grouping value of the data, which allows overlain networks, with the default being 'pol'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @param dest This is the destination folder where the function output will be saved. Defaults to current working directory.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10'.
#' @param long.term This sets the length of the long-term trend. Defaults to 1 week.
#' @param short.term This sets the length of the short-term trend. Defaults to 1 day.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' drFUN()

drFUN <- function(x, obs, group1, group2, reflective = TRUE, plot = TRUE, dest = getwd(), date.start = '2016-07-01', date.end = '2016-07-10', long.term = 168, short.term = 24, theta = NA, tau = NA) {
 
	# install and load required packages
		list.of.packages <- c("stats", "stringr", "raster", "ggplot2", "data.table", "dplyr")
		new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
		if(length(new.packages)) install.packages(new.packages)
		library(stats); library(stringr); library(raster); library(ggplot2); library(data.table); library(dplyr)
		
	# clause on type of analysis to be run, define dates
	if (reflective == TRUE) {
	
		date.start = date.start 
		date.end = date.end
	} else {
	
		date.start = str_c(Sys.Date()-7)
		date.end = str_c(Sys.Date())
	}
	
	# define selected variables
	# use data.table package to deal with large datasets
		setDT(x)
		x[, date := ymd(date)]
		x <- x[date >= date.start & date <= date.end]
		x$group1 <- x[,..group1]
		x$group2 <- x[,..group2]
		x$obs <- x[,..obs]
	
	# training data - select first [long.term] length of data
		z.train <- x %>%
			head(n = long.term)
  
	# derive the average
		z.train.av <- mean(z.train[, dat], na.rm = T)
  
	# remove the average from the observations & force NA to 0 for matrix calculations
		x.obs <- z.train[, dat] - z.train.av
		x.obs[is.na(x.obs)] <- 0
  
	# short-term autocorrelation matrix to [short.term] length
		x.obs.lag <- embed(x.obs, short.term)
  
	# X here is the data to be predicted, H here is the data used for predictions
		x.obs.X <- x.obs.lag[, 1]
		x.obs.H <- x.obs.lag[, short.term:2]
  
	# matrix coefficients (N-1 length)
		b <- solve(qr(x.obs.H), x.obs.X)
   
	# application of above to the following data
	# derive M length averages - these will be rolling at each unit
	# the first M data will use the set average
  
		df <- x %>%
			mutate(movingav = movingFun(x[, dat], long.term, fun = mean, type = 'to', na.rm = T)) %>%
			mutate(movingav = replace(movingav, 1:long.term, z.train.av)) %>%
			mutate(X = x[, dat] - movingav)
  
	# set up the matrix
		x.obs.lag <- embed(df$X, short.term)
		x.obs.X <- x.obs.lag[, 1]
		x.obs.H <- x.obs.lag[, short.term:2]
  
	# solve the matrix using the derived coefficients from the model
		b.X <- sweep(x.obs.H, MARGIN = 2, b, '*')
  
	# this sums all coefficient predictions at each unit to give an estimate of column N+1
		x.Hb <- rowSums(b.X, na.rm = T)

	# determine the difference between observed and predicted
		iX <- x.obs.X - x.Hb

	# insert some leading N values
		rep.NA <- rep(NA, short.term-1)
		iX <- c(rep.NA, iX)
	
	# join results to original data file
		iX.data <- cbind(x, iX) 
     
	# use theta and tau thresholds if present
	if((!is.na(theta))) {
  
		if((!is.na(tau))) {
  
	# if else clause on whether the data are greater than theta, becomes 1 if true, and 0 otherwise
	# this looks at running means and is for tau, the length of time for alarms
	
			iX.data <- iX.data[, change.thresh := ifelse(iX > theta, 1, 0), by = list(group1, group2)][, alarm := movingFun(change.thresh, n = tau, type = 'to', fun = mean, na.rm = T)] 
			 
		} else {
		
		tau <- NA
		} 
	}  else {
	
	theta <- NA
	tau <- NA
	}
  
	# create csv file to be saved to dest folder (default is directory)
		csv <- write.csv(iX.data, file = str_c(dest, '/decomposedRegression', start.date, '.csv'))
  
	# create plot and save
		if(plot == TRUE) {
  
	# create plot of the ranked variance over the last week
		plot.name <- str_c(dest, '/decomposedRegression', start.date, '.pdf')
		pdf(plot.name, width = 6, height = 10)
		ts <- suppressWarnings(ggplot(iX.data, aes(x = as.POSIXct(day), y = iX, colour = group2)) +
					theme_bw() +
					geom_point() +
					facet_wrap(~ group1) +
					xlab("") + ylab('residual'))
  
		if(!is.na(theta)) {ts <- ts + geom_abline(intercept = theta, slope = 0)}
	dev.off()
	
	return(c(ts, csv))
	} else { 
  
	# returns x as saved csv file
		csv
  }
}


