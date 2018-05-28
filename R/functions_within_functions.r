#' A selection of functions used within other functions in the package
#' 
#' A function for creating dissimilarity plots using Spearman rank correlations.
#' @param omit.cols This is a list of the columns to remove from analysis if needed. Default is 
#' @export
#' @examples
#' distFUN()

distFUN <- function(x, omit.cols = NA) {
  
  library(stringr);library(lubridate);
  
  x <- as.data.frame(x);

  # clauses for where columns are empty, all zeros or same number, etc.
  row.names <- colnames(x)[!(colnames(x) %in% omit.cols)]; 
  xx <- sapply(x[, row.names], var, na.rm = T) == 0;
  xx <- which(xx == T);
    
  if(length(!is.na(xx))){
    xxx <- x[, -which(names(x) %in% names(xx))];
    row <-row.names[! row.names %in% names(xx)];
  } else {
    xxx <- x;
    row <- row.names;
  };
    
  # identify empty observations
  yy <- colnames(x[, row.names])[colSums(x[, row.names], na.rm = T) == 0];   
  if(length(yy > 0)){      
    xxx <- xxx[, -which(names(xxx) %in% yy)];
    row <- row.names[! row.names %in% yy];	
  };
    
  # apply dissimilarity correlation matrix
  row.names <- colnames(xxx)[!(colnames(xxx) %in% omit.cols)];
  dist <- as.matrix(1 - abs(cor(xxx[, row.names], use = 'pairwise.complete.obs', method = 'spearman')));
    
  # this is to identify columns which contain NA, etc - need to remove specific moments
  zz <- names(which(colSums(is.na(dist)) > 0.75 * length(rownames(dist))));
    
  # removal of those observations where over 75% of comparisons are NA
  if(length(zz > 0)){
    dist <- dist[-which(row.names(dist) %in% zz), -which(row.names(dist) %in% zz)];
    row <- row.names[! row.names %in% zz];
  }
  
  # turn NA results into 0 - not ideal, however appears to be best option
  dist[is.na(dist)] <- 0;
    
  # convert back to dissimilarity structure
  # apply a hierarchical cluster function
  hc <- hclust(as.dist(dist));
  
  # this finds the dissimilarity score and where an instrument first joins a cluster
  hc.labels <- data.frame(list = 1:length(hc$labels), labels = hc$labels);
  hc.times <- data.frame(order = 1:length(hc$merge[,1]), from = ifelse(hc$merge[,1] < 0, -1*hc$merge[,1], 0), to = ifelse(hc$merge[,2] < 0, -1*hc$merge[,2], 0));
  hc.times <- melt(hc.times, id.vars = 'order');
  hc.times <- hc.times %>% filter(value != 0);
  hc.times$list <- hc.times$value; hc.times$variable <- NULL; hc.times$value <- NULL;
  hc.labels.time <- join(hc.labels, hc.times, by = 'list');
  hc.labels.time <- hc.labels.time %>% arrange(order);
   
  hc.from <- ifelse(hc$merge[,1] < 0, -1, 1);
  hc.to <- ifelse(hc$merge[,2] < 0, -1, 1);
  hc.sum <- hc.from + hc.to;
  hc.sum[hc.sum == 0] <- 1;
  hc.sum[hc.sum == 2] <- 0;
  hc.sum[hc.sum == -2] <- 2;

  final.height <- rep(hc$height, hc.sum);
  hc.labels.time$height <- final.height;
  
  f <- hc.labels.time %>%
	dplyr::select(group = labels, dissimilarity = height);
	
};
  

#' This function gives the proxy type of network median.
#'
#' @param obs This is the column being examined.
#' @param group This is the identifying variable.
#' @param id This is an option if some identifying column is required (e.g. by site). Default is NA.
#' @export
#' @examples
#' networkmedianFUN()

networkmedianFUN <- function(x, group, obs, by.row = T, id = NA, statistic = median){
  
  library(data.table); library(dplyr);

  # define variables
  x <- as.data.table(x);
  x$obs <- x[[obs]]; x$group <- x[[group]];
  if(!is.na(id)) x$id <- x[[id]];
  
  # filter data to that of interest
  z = x[, list(group, obs)];
  if(!is.na(id)) z.id <- x[, list(group, id, obs)];
  
  net.day.FUN <- function(z){
    if(length(z$obs) > 1){
      network.proxy <- unlist(sapply(seq(1, nrow(z)), function(i){
        test = z[-i,]; 
        network.proxy <- statistic(test$obs, na.rm = T);
        return(network.proxy);
      }))
      } else {
        network.proxy <- NA;
        network.proxy <- as.numeric(network.proxy);
      }
    z <- cbind(z, network.proxy);
    z <- unique(z);
    };
  
  net.FUN <- function(z){
    group.list <- as.vector(unique(z$group));
    output = data.frame(group = as.character(),
                        network.proxy = as.numeric());
    for(i in group.list){
    if(length(z$obs) > 1){
      selected.group = i;
      group.stat <- data.frame(group = selected.group, 
                               network.proxy = statistic(z$obs[z$group != selected.group], na.rm = T));
      output = rbind(output, group.stat);
    }
    }
    x <- full_join(x, output);
    x <- unique(x);
  }
  
  if(by.row == T){
    Z.data <- z[, net.day.FUN(.SD), by = group];
  }
  
  if(by.row == F){
    Z.data <- net.FUN(z)
  }
  
  if(!is.na(id)) Z.data <- full_join(Z.data, z.id);
  return(Z.data);
};

#' This function gives the proxy type of nearest site data.
#'
#' @param obs This is the column being examined.
#' @param group This is the grouping variable.
#' @param lat Latitude column.
#' @param lon Longitude column.
#' @export
#' @examples
#' nearestsiteFUN()
  
nearestsiteFUN <- function(x, obs, group, lat, lon){

  library(sp);library(rgeos);library(dplyr);

  # define variables
  x <- as.data.frame(x);
  x$obs <- x[, obs]; x$group <- x[, group]; x$lat <- x[, lat]; x$lon <- x[, lon];
  
  # find closest sites
  y <- x[group %in% unique()];
    meta %>%
    filter(group %in% unique(x$group)) %>%
    select(group, lat, lon);
  
  sp.y <- y;	
  coordinates(sp.y) <- ~lon+lat;
  distances <- gDistance(sp.y, byid = T);	
  min.distance <- apply(distances, 1, function(x) order(x, decreasing = F)[2]);
  nearest.dist <- cbind(y, y[min.distance,], apply(distances, 1, function(x) sort(x, decreasing=F)[2]));
  colnames(nearest.dist) <- c(colnames(y), 'n.group', 'n.lat', 'n.lon', 'distance');

  # combine the data with appropriate proxy
  df3 <- full_join(x, nearest.dist, by = 'group');
  x$nearest.proxy <- x$obs;
  x$obs <- NULL;
	
  finaldata <- left_join(df3, x, by = c('date', 'n.group' = 'group'));
  finaldata <- finaldata %>%
    select(date, group, obs, nearest.proxy);
	
}


#' This function gives the rolling KS function
#'
#' @param obs This is the column being examined.
#' @param proxy This is the comparison variable
#' @param window This defines the size of the window to sample data from 
#' @export
#' @examples
#' rollingKStest()

rollingKStest <- function(z, obs, proxy, window = 1440){
  
  library(zoo);library(dplyr);

  # define variables
  z <- as.data.frame(z);
  z$obs <- z[, obs]; z$proxy <- z[, proxy];
  
  if(length(z$obs) > window){
    
    z <- z %>%
      dplyr::arrange(date);
    
    # save date column
    date <- z[, 'date'];
    
    # turn data into data frame
    y <- data.frame(z[, obs], z[,proxy]);
  
    # convert df into zoo classes used in R for rolling functions
    y.zoo <-zoo(y);
    colnames(y.zoo) <- c('r', 's');
  
	# create custom KS function to deal with missing data
	
	ks.dev.p <- function(x, y){
  
  min.length <- min(c(length(na.omit(x)),length(na.omit(y))));
  
  if(min.length > 0.5 * length(x)){
    p.value <- ks.test(x, y)$p.value;
  } else {
    p.value <- NA;
  }
  return(p.value);
};

ks.dev.s <- function(x, y){
  
  min.length <- min(c(length(na.omit(x)),length(na.omit(y))));
  
  if(min.length > 0.5 * length(x)){
    statistic <- ks.test(x, y)$statistic;
  } else {
    statistic <- NA;
  };
  return(statistic);
};
    # running ks functions
    p.value <- rollapply(y.zoo, function (x) ks.dev.p(x[, 'r'], x[, 's']), 
	  by.column = F, fill = NA, align = 'right', width = window);
    statistic <- rollapply(y.zoo, function (x) ks.dev.s(x[, 'r'], x[, 's']), 
	  by.column = F, fill = NA, align = 'right', width = window);
  
    p.value <- data.frame(p.value);
    statistic <- data.frame(statistic);
  
    # make results into df and tidy
    model <- cbind(p.value, statistic);
    names(model) <- c('p.value', 'statistic');
    model <- data.frame(lapply(model, function(x) round(x, 6)));
    model <- cbind(date, model);
  
    # join ks results to the data and return
    z <- join(z, model, by = 'date');
    } else {	
	
	z$p.value <- NA;
	z$statistic <- NA;
 };
};
