#' A selection of functions used within other functions in the package
#' 
#' A function for creating dissimilarity plots using Spearman rank correlations.
#' @param date Date column
#' @export
#' @examples
#' distFUN()

distFUN <- function(x, date) {
  
  # install and load required functions  
  list.of.packages <- c("stats", "stringr", "lubridate")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(stats); library(stringr); library(lubridate)
  
  x <- as.data.frame(x)
  
  # find date
	date <- x[, date]
  
  # clauses for where columns are empty, all zeros or same number, etc.
  row.names <- colnames(x[, -c(1, ncol(x))])
      
  xx <- sapply(x[, row.names], var, na.rm = T) == 0
  xx <- which(xx == T)
    
  if(length(!is.na(xx))){
    xxx <- x[, -which(names(x) %in% names(xx))]
    row <-row.names[! row.names %in% names(xx)]
  } else {
    xxx <- x
    row <- row.names
  }
    
	# identify empty observations
  yy <- colnames(x[, row.names])[colSums(x[, row.names], na.rm = T) == 0]
    
  if(length(yy > 0)) {      
    xxx <- xxx[, -which(names(xxx) %in% yy)]
    row <- row.names[! row.names %in% yy]	
  }
  
  # if clause on days with less than values 

  
  # apply dissimilarity correlation matrix
  dist <- as.matrix((1 - cor(xxx[, -c(1, ncol(xxx))], use = 'pairwise.complete.obs', 
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
  
  hc.labels <- data.frame(list = 1:length(hc$labels), labels = hc$labels)
  hc.times <- data.frame(order = 1:length(hc$merge[,1]), from = ifelse(hc$merge[,1] < 0, -1*hc$merge[,1], 0), to = ifelse(hc$merge[,2] < 0, -1*hc$merge[,2], 0))
  hc.times <- melt(hc.times, id.vars = 'order')
  hc.times <- hc.times %>% filter(value != 0)
  hc.times$list <- hc.times$value; hc.times$variable <- NULL; hc.times$value <- NULL
  hc.labels.time <- join(hc.labels, hc.times, by = 'list')
  hc.labels.time <- hc.labels.time %>% arrange(order)
   
  hc.from <- ifelse(hc$merge[,1] < 0, -1, 1)
  hc.to <- ifelse(hc$merge[,2] < 0, -1, 1)
  hc.sum <- hc.from + hc.to
  hc.sum[hc.sum == 0] <- 1
  hc.sum[hc.sum == 2] <- 0
  hc.sum[hc.sum == -2] <- 2

  final.height <- rep(hc$height, hc.sum)

  hc.labels.time$height <- final.height
  
  f <- hc.labels.time %>%
	select(site = labels, dissimilarity = height)
	
}
  

#' This function gives the proxy type of network median.
#'
#' @param obs This is the column being examined.
#' @param group This is the identifying variable.
#' @export
#' @examples
#' networkmedianFUN()

networkmedianFUN <- function(x, group, obs){
  
  list.of.packages <- c("stats", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(stats); library(dplyr)

  # define variables
  x <- as.data.table(x)
  x$obs <- x[, ..obs]; x$group <- x[, ..group]
  
  # filter data to that of interest
  z = x[, list(group, obs)]
  
  if(length(z$obs) > 1){
    network.proxy <- unlist(sapply(seq(1, nrow(z)), function(i) {
      test = z[-i,]  
      z <- median(test$obs, na.rm = T)  
      return(z)
  }))
  } else {
    network.proxy <- NA
    network.proxy <- as.numeric(network.proxy)
  }
  
  z <- cbind(z, network.proxy)
  z <- unique(z)
}

#' This function gives the proxy type of nearest site data.
#'
#' @param obs This is the column being examined.
#' @param group This is the grouping variable.
#' @param lat Latitude column.
#' @param lon Longitude column.
#' @export
#' @examples
#' nearestsiteFUN()
  
nearestsiteFUN <- function(x, obs, group, lat, lon) {

  # install and load required packages
  list.of.packages <- c("dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(dplyr)

  # define variables
  x <- as.data.frame(x)
  x$obs <- x[, obs]; x$group <- x[, group]; x$lat <- x[, lat]; x$lon <- x[, lon]
  
  # find closest sites
  y <- x[group %in% unique()]
    meta %>%
    filter(group %in% unique(x$group)) %>%
    select(group, lat, lon)
	  
  d <- as.matrix(dist(cbind(x$lon, x$lat)))
  d <- data.frame(d)
  colnames(d)<- c(as.character(d$group), 'group')	
  min.d <- apply(d, 1, function(x) order(x, decreasing = F)[2])
  newdata <- cbind(meta2, meta2[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))
  colnames(newdata) <- c(colnames(meta2), 'n.group', 'n.lat', 'n.lon', 'distance')
  newdata <- newdata[,c(1, 5:9)]

  # combine the data with appropriate proxy
  df3 <- join(x, newdata, by = 'group')
  x$nearest.proxy <- x$obs
  x$obs <- NULL
	
  finaldata <- left_join(df3, x, by = c('date', 'n.group' = 'group'))
  finaldata <- finaldata %>%
    select(date, group, obs, nearest.proxy)
	
}


#' This function gives the rolling KS function
#'
#' @param obs This is the column being examined.
#' @param proxy This is the comparison variable
#' @param window This defines the size of the window to sample data from 
#' @export
#' @examples
#' rollingKStest()

rollingKStest <- function(z, obs, proxy, window = 1440) {
  
  list.of.packages <- c("zoo", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  library(zoo); library(dplyr)

  # define variables
  z <- as.data.frame(z)
  z$obs <- z[, obs]; z$proxy <- z[, proxy]
  
  if(length(z$obs) > window){
    
    z <- z %>%
      arrange(date)
    
    # save date column
    date <- z[, 'date']
    
    # turn data into data frame
    y <- data.frame(z[, obs], z[,proxy])
  
    # convert df into zoo classes used in R for rolling functions
    y.zoo <-zoo(y)
    colnames(y.zoo) <- c('r', 's')
  
	# create custom KS function to deal with missing data
	
	ks.dev.p <- function(x, y){
  
  min.length <- min(c(length(na.omit(x)),length(na.omit(y))))
  
  if(min.length > 0.5 * length(x)){
    p.value <- ks.test(x, y)$p.value
  } else {
    p.value <- NA
  }
  return(p.value)
}

ks.dev.s <- function(x, y){
  
  min.length <- min(c(length(na.omit(x)),length(na.omit(y))))
  
  if(min.length > 0.5 * length(x)){
    statistic <- ks.test(x, y)$statistic
  } else {
    statistic <- NA
  }
  return(statistic)
}
    # running ks functions
    p.value <- rollapply(y.zoo, function (x) ks.dev.p(x[, 'r'], x[, 's']), 
	  by.column = F, fill = NA, align = 'right', width = window)
    statistic <- rollapply(y.zoo, function (x) ks.dev.s(x[, 'r'], x[, 's']), 
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