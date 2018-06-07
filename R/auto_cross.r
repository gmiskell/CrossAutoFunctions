#' A function to set up data ready to perform "auto" analysis.
#'
#' This function reformats the data relative to auto analysis, where site is held constant.
#' @param obs This is the column under observation.
#' @param date This is the data column.
#' @param group This is the id of each set of data being compared.
#' @param ell This sets the length of time backwards (in days) that should be considered as the auto proxy. Defaults to 3.
#' @param sample.unit This classifies whether the analysis should consider data lags by days ('day', default), or by hours ('hour').
#' @export
#' @examples
#' autoFUN()

autoFUN <- function(x, date, obs, ell = 3, sample.unit = 'day', group = 'site'){
  
library(plyr);library(raster);library(data.table);library(stringr);library(lubridate);library(dplyr);
  
	x <- as.data.table(x);
	x$group <- x[[group]]; x$obs <- x[[obs]]; x$date <- x[[date]];
	
	x <- x[, list(date, group, obs)];
	x <- unique(x);
	x <- x[order(group)];
	
	len <- as.numeric(difftime(x$date[2], x$date[1], units = 'mins'));
	if(len == 60) len = 24;
	if(len == 1) len = 1440;
	
	# pad out of data frame in case of missing values
	if(sample.unit == 'day') {
	  if(len == 24) time.freq = 'hour' 
	  if(len == 1440) time.freq = 'min'
	  nn <- seq(from = 0, to = ell) * len
	};
	
	if(sample.unit == 'hour'){
	  if(len == 24) {ell = len * ell; time.freq = 'hour'}
	  if(len == 1440) {ell = (len * ell) / 60; time.freq = 'min'}
	  nn <- seq(from = 0, to = ell)	
	};
	
	span <- as.POSIXct(seq(floor_date(min(x$date, na.rm = T), unit = 'day'), max(x$date, na.rm = T), by = time.freq));
	grp.n <- unique(na.omit(x$group));
	span <- rep(span, length(grp.n));
	grp <- rep(grp.n, each = length(span));
	span <- data.frame(date = span, group = grp);
	
	x <- merge(span, x, by = c('group', 'date'), all = TRUE);
	x <- unique(x);
  
	x <- setDT(x)[, paste("l", 0:ell) := shift(obs, n =  nn, type = 'lag'), by = group];

	y <- x[, melt(x, id.vars = c('date', 'group', 'obs'))][!is.na(group)];
	setnames(y, c('variable', 'value'), c('comparison', 'comparison.value'));
  
	y;

}

#' A function to set up data ready to perform "cross" analysis.
#'
#' This function reformats the data relative to cross analysis, where time is held steady. 
#' @param obs This is the column under observation.
#' @param group This is the id of each set of data being compared.
#' @export
#' @examples
#' crossFUN()

crossFUN <- function(x, date, obs, group){
  
	library(reshape2);library(lubridate);library(data.table);
  
	x <- as.data.table(x);
	
	x$group <- x[[group]]; x$obs <- x[[obs]]; x$date <- x[[date]];
	x$group <- as.character(x$group); x$site <- as.character(x$site);
	
	x <- x[, list(date, group, obs)][, date := ymd_hms(date)];
	cast.x <- x[, dcast(x, ... ~ group, median, value.var = 'obs')];

	y <- join(x, cast.x, by = 'date');
	y <- y[, melt(y, id.vars = c('date', 'group', 'obs'))];
	setnames(y, c('variable', 'value'), c('comparison', 'comparison.value'));
	y <- y[!(y$group == y$comparison),];
	
	y;
	
}
