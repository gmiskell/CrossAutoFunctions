#' A selection of functions used within other functions in the package
#' 
#' A function for creating dissimilarity plots using Spearman rank correlations.
#' @param d.thresh This is the dissimilarity threshold (1-p).
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @export
#' @examples
#' distFUN()

distFUN <- function(x, d.thresh = 0.3, plot = TRUE) {
  
  # install and load required functions  
  list.of.packages <- c("stats", "stringr", "lubridate")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(stats); library(stringr); library(lubridate)
  
  # clauses for where columns are empty, all zeros or same number, etc.
  row.names <- colnames(x[, -c(1:2, ncol(x))])
  
  x$day <- day(x$date)
    
  xx <- sapply(x[, row.names], var, na.rm = T) == 0
  xx <- data.frame(which(xx))
    
  if(length(!is.na(xx))){
    xxx <- x[, -which(names(x) %in% rownames(xx))]
    row <-row.names[! row.names %in% rownames(xx)]
  } else {
    xxx <- x
    row <- row.names
  }
    
  yy <- colnames(x[, row.names])[colSums(x[, row.names], na.rm = T) == 0]
    
  if(length(yy > 0)) {      
    xxx <- xxx[, -which(names(xxx) %in% yy)]
    row <- row.names[! row.names %in% yy]	
  }
    
  # apply dissimilarity correlation matrix
  dist <- as.matrix((1 - cor(xxx[, -c(1:2, ncol(xxx))], use = 'pairwise.complete.obs', 
                           method = 'spearman'))/ 2)
    
  # this is to identify columns which contain NA, etc - need to remove specific moments
  zz <- names(which(colSums(is.na(dist))>0.75*length(rownames(dist))))
    
  # removal of those observations where over 75% of comparisons are NA
  if(length(zz > 0)){
    dist <- dist[-which(row.names(dist) %in% zz), -which(row.names(dist) %in% zz)]
    row <- row.names[! row.names %in% zz]
  }
  
  # turn NA results into 0 - not ideal, however appears to be best option
  dist[is.na(dist)] <- 0
    
  # convert back to dissimilarity structure
  # apply a hierarchical cluster function
  hc <- hclust(as.dist(dist))
  
	
  # create final output
  if(plot == TRUE){ 
	plot(hc, xlab = str_c(x$day[1], x$type[1]), ylab = 'dissimilarity')
	abline(a=d.thresh, b=0)
	} else {
	return(dist)
	}
	
}
  

#' This function gives the proxy type of network median
#' @param obs This is the column being examined.
#' @param site This is the grouping variable
#' @param type A second grouping variable
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @export
#' @examples
#' networkmedianFUN()

networkmedianFUN <- function(x, obs, site, type){
  
  list.of.packages <- c("stats", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(stats); library(dplyr)

  # filter data to that of interest, removing empty rows
  z = x %>%
    select(date, obs, site, type) %>%
    na.omit()
  
  if(length(z[[1]]) < 2000) {
  # make comparison data
  proxy <- unlist(sapply(seq(1, nrow(z)), function(i) {
    test = z[-i,]  
    proxy = median(test$obs, na.rm = T)   
  }))} else {
    proxy = z[, median(obs, na.rm = T)]
  }
}


#' This function gives the proxy type of nearest site data
#' @param obs This is the column being examined.
#' @param site This is the grouping variable
#' @param type A second grouping variable
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#' @export
#' @examples
#' nearestsiteFUN()
  
nearestsiteFUN <- function(x, meta, site, type) {

  # install and load required packages
  list.of.packages <- c("dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(dplyr)

  # label current pol
  pol.lab <- x$type[1]
	
  # find closest sites
  meta2 <- meta %>%
    filter(type == pol.lab) %>%
    filter(site %in% unique(x$site)) %>%
    select(site, type, lat, lon)
	  
  d <- as.matrix(dist(cbind(meta2$lon, meta2$lat)))
  d <- data.frame(d)
  d$site <- meta2$site
  colnames(d)<- c(as.character(d$site), 'site')	
  min.d <- apply(d, 1, function(x) order(x, decreasing = F)[2])
  newdata <- cbind(meta2, meta2[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))
  colnames(newdata) <- c(colnames(meta2), 'n.site', 'n.pol', 'n.lat', 'n.lon', 'distance')
  newdata <- newdata[,c(1, 5:9)]

  # combine the data with appropriate proxy
  df3 <- join(x, newdata, by = 'site')
  x$proxy <- x$obs
  x$obs <- NULL
	
  finaldata <- left_join(df3, x, by = c('date', 'n.site' = 'site'))
  finaldata <- finaldata %>%
    select(site, date, type = type.x, obs, proxy)
	
}


#' This function gives the rolling KS function
#' @param obs This is the column being examined.
#' @param proxy This is the comparison variable
#' @param window. This defines the size of the window to sample data from 
#' @export
#' @examples
#' rollingKStest()

rollingKStest <- function(z, obs, proxy, window = 1440) {
  
  list.of.packages <- c("zoo")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(zoo)

  if(length(z[, obs] > window)){
    
    z <- z %>%
      arrange(date)
    
    # save date column
    date <- z[, 'date']
    
    # turn data into dataframe
    y <- data.frame(z[, obs], z[,proxy])
  
    # convert df into zoo classes used in R for rolling functions
    y.zoo <-zoo(y)
    colnames(y.zoo) <- c('r', 's')
  
    # running ks functions
    p.value <- rollapply(y.zoo, function (x) ks.test(x[, 'r'], x[, 's'])$p.value, 
	  by.column = F, fill = NA, align = 'right', width = window)
    statistic <- rollapply(y.zoo, function (x) ks.test(x[, 'r'], x[, 's'])$statistic, 
	  by.column = F, fill = NA, align = 'right', width = window)
  
    p.value <- data.frame(p.value)
    statistic <- data.frame(statistic)
  
    # make results into df and tidy
    model <- cbind(p.value, statistic)
    names(model) <- c('p.value', 'statistic')
    model <- data.frame(lapply(model, function(x) round(x, 6)))
    model <- cbind(date, model)
  
    # join ks results to the data and return
    z <- join(z, model, by = 'date')
 } else {	
	z$p.value <- NA
	z$statistic <- NA
 }
  
}


#' This function provides the KS plot.
#' @param obs This is the column being examined.
#' @param proxy This is the comparison variable
#' @param site This is the grouping factor
#' @export
#' @examples
#' ksplotFUN()

ksplotFUN <- function(x, obs, proxy, site){
  	
  list.of.packages <- c("ggplot2", "gridExtra")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
 
  library(ggplot2); library(gridExtra)

  label.1min <- ks.test(x$obs[x$time == '1 min'], x$proxy[x$time == '1 min'])$p.value
  label.1hr <- ks.test(x$obs[x$time == '1 hr'], x$proxy[x$time == '1 hr'])$p.value
    
  labels = data.frame(time = c('1 min', '1 hr'), ks = c(round(label.1min, 3), round(label.1hr, 3)))
	
  plotA <- ggplot(x, aes(x = obs, colour = 'obs.')) + theme_bw() +
    stat_ecdf() + stat_ecdf(aes(x = proxy, colour = 'proxy')) + facet_wrap(~ time, scales = 'free_x') + 
	xlab("ranked concentration") + ylab("Fn(X)") + theme(legend.title = element_blank(), legend.position = 'bottom') +
	geom_text(data = labels, aes(label = paste("KS test: ", ks), x = 5, y = 0.15, parse = T), colour = 'black')
  plotB <- ggplot(x, aes(x = date, y = p.value)) + theme_bw() + geom_point() + facet_wrap(~ time) + ylim(c(0,1))
	
  grid.arrange(plotA, plotB, top = x$site[1])
	
}