#' A decomposed regression analysis.
#'
#' This function gives the decomposed time series with long-term (week) and short-term (day) trends removed.
#' @param obs This is the numeric column under investigation.
#' @param proxy This is the data that obs is being compared to (cannot be the same as obs).
#' @param date This is a required column with the format 'YYYY-MM-DD HH:MM:SS'.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param date.start This sets the start of the selected data, if reflective is true. Defaults to '2016-07-01 00:00:00'.
#' @param date.end This sets the end of the selected data, if reflective is true. Defaults to '2016-07-10 00:00:00'.
#' @param long.term This sets the length of the long-term trend. Defaults to 168.
#' @param short.term This sets the length of the short-term trend. Defaults to 24.
#' @param theta This sets the test threshold on whether to flag the data. Defaults to NA (no flags given).
#' @param tau This sets the day threshold on whether to flag the data given consistent theta flags. Defaults to NA (no flags given).
#' @export
#' @examples
#' drFUN()

drFUN <- function(x, obs, date, proxy, reflective = TRUE, date.start = '2016-07-01 00:00:00', date.end = '2016-07-10 00:00:00', long.term = 168, short.term = 24, theta = NA, tau = NA) {
 
	# install and load required packages
		list.of.packages <- c("stats", "stringr", "raster", "ggplot2", "data.table", "dplyr")
		new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
		if(length(new.packages)) install.packages(new.packages)
		library(stats); library(stringr); library(raster); library(ggplot2); library(data.table); library(dplyr)
		
		# clause on type of analysis to be run
		if (reflective == TRUE) {
		  
		  date.start = date.start 
		  date.end = date.end
		} else {
		  
		  date.start = Sys.time()-60*60*24*7
		  date.end = Sys.time()
		}
		
		date.interval <- interval(date.start, date.end, tzone = 'Pacific/Auckland')

		# define selected variables
		# use `data.table` package to deal with large datasets
		x <- as.data.table(x)
		x[, date := ymd_hms(date, tz = 'Pacific/Auckland')]
		x <- x[date %within% date.interval][, date := force_tz(date, tzone = 'Pacific/Auckland')]
		x$proxy <- x[,..proxy]
		x$obs <- x[,..obs]
	
		# check length of data
		if(length(na.exclude(x$obs)) > 0.8*long.term){	
		
	# training data - select first [long.term] length of data using the proxy data
		z.train <- head(x, n = long.term)
  
	# derive the average
		z.train.av <- mean(z.train[, proxy], na.rm = T)
  
	# remove the average from the observations & force NA to 0 for matrix calculations
		x.obs <- z.train[, proxy] - z.train.av
		x.obs[is.na(x.obs)] <- 0
  
	# short-term autocorrelation matrix to [short.term] length
		x.obs.lag <- embed(x.obs, short.term)
  
	# X here is the data to be predicted, H here is the data used for predictions
		x.obs.X <- x.obs.lag[, 1]
		x.obs.H <- x.obs.lag[, short.term:2]
  
	# matrix coefficients (N-1 length)
		b <- solve(qr(x.obs.H), x.obs.X)
   
	# application of above to the obs (Y) and proxy (X) data
	# derive [long.term] length averages
	# the first [long.term] data will use the model average
		uni.date <- x[,unique(date)]
		breaks <- gl(ceiling((length(uni.date))), long.term)[1:(length(uni.date))]
		date.breaks <- data.table(date = uni.date, breaks)
		
		df <- join(x, date.breaks, by = 'date')
		df <- setDT(df)[, lag := lag(proxy, long.term)]
		df <- df[, lt.av := mean(lag, na.rm = T), by = breaks][, lt.av := ifelse(lt.av == 'NaN', z.train.av, lt.av)][, X := obs - lt.av][, Y := proxy - lt.av]
  
	# set up the matrix
		prxy.lag <- embed(df$X, short.term)
		prxy.X <- prxy.lag[, 1]
		prxy.H <- prxy.lag[, short.term:2]
		obs.lag <- embed(df$Y, short.term)
		obs.Y <- obs.lag[, 1]
		obs.H <- obs.lag[, short.term:2]
  
	# solve the matrix using the derived coefficients from the model
		b.X <- sweep(obs.H, MARGIN = 2, b, '*')
		b.Y <- sweep(prxy.H, MARGIN = 2, b, '*')
  
	# this sums all coefficient predictions at each unit to give an estimate of column N+1
		x.Hb <- rowSums(b.X, na.rm = T)
		y.Hb <- rowSums(b.Y, na.rm = T)

	# determine the difference between observed and predicted
		iX <- prxy.X - x.Hb
		iY <- obs.Y - y.Hb

	# insert some leading N values
		rep.NA <- rep(NA, short.term-1)
		iX <- c(rep.NA, iX)
		iY <- c(rep.NA, iY)
		xy.diff <- iX - iY
	
	# join results to original data file
		XY.data <- cbind(x, iX, iY, xy.diff) 
     
	# use theta and tau thresholds if present
	  if((!is.na(theta))) {
  
	    XY.data <- setDT(XY.data)[, warning := ifelse(abs(xy.diff) > theta, 1, 0)]
	  
	  	if((!is.na(tau))) {
  
	# if else clause on whether the data are greater than theta, becomes 1 if true, and 0 otherwise
	# this looks at running means and is for tau, the length of time for alarms
	
		    XY.data <- XY.data[, alarm := movingFun(warning, n = tau, type = 'to', fun = mean, na.rm = T)] 
		  } else {
		
		    XY.data$alarm <- NA
		  } 
	  }  else {
	
	  XY.data$warning <- NA
	  XY.data$alarm <- NA
	  }

	# gather variables of interest and return
		XY.data <- XY.data[, list(date, test = 'decomposed regression', statistic = xy.diff, warning, alarm)]
		} else {
		  
		  XY.data <- x[, list(date, test = 'decompose regression', statistic = NA, warning = NA, alarm = NA)]
		  }

	return(XY.data)
}


