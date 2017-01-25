#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' idwcFUN()

# inverse distance correlation filter script

# description:  this filter will compare the correlations of a site against the other 
#               sites. a weighting will be applied, giving higher value to those sites 
#               closer. the final product will be a sum of the correlations multiplied 
#               by the inverse distance

idwcFUN <- function(df, meta = meta.data, correlation.threshold = NA, 
                                  window.length = 1440){
  
    # load libraries
    library(openair);library(plyr);library(stringr);
    library(reshape2);library(dplyr)


    # calculate distances among locations -- this could be saved as a table prior
    distances <- as.matrix(dist(cbind(meta$lon, meta$lat)))
    distance <- data.frame(distances)
    distance$site <- str_c(meta$site, '_', meta$pol)
    # remove empty locations or any errors in the data
    distance <- distance[!(is.na(distance$site) | distance$site==""), ]
    distance <- Filter(function(x)!all(is.na(x)), distance)
    colnames(distance)<- c(as.character(distance$site), 'site')
    # rearrange data frame and convert results to km
    dist.melt <- distance %>%
      melt(id.vars = 'site') %>%
      mutate(distance = value*100) %>%
      mutate(variable = str_c(site, '.', variable)) %>%
	  select(-value, -site) 

	
    # rearrange the data to perform correlations per column
    df.3 <- df %>%
      select(date, site, pol, conc_mean) %>%
      dcast(... ~ site + pol, median)
	

    # create sequence with blocks of length window to perform analysis upon
	breaks <- rep(1:ceiling(nrow(df.3)/window.length), each = window.length)
	breaks <- breaks[1:nrow(df.3)]
	
		
	# perform correlation on the data for block size
	break.cor <- ddply(df.3, .(breaks),
	  function(x) cor(x[,-1], use = 'pairwise.complete.obs', method = 'spearman'))
	names <- colnames(break.cor[,-1])
	break.cor$site <- rep(names, length(unique(break.cor$breaks)))
	break.cor$site2 <- rep(names, each = length(unique(break.cor$breaks)))
	

    # join distance and correlations and create filter
    corIDW <- break.cor %>%
      melt(id.vars = c('breaks', 'site', 'site2')) %>%
	  mutate(variable = str_c(site, '.', variable)) %>%
	  select(-site, -site2) %>%
      rename(correlation = value) %>%
      join(dist.melt, by = 'variable') %>%
      # remove correlations performed on the same site and NA correlations
      filter(distance > 0 & !is.na(correlation)) %>%
      # add in the inverse distance
      mutate(filt = correlation * 1 / distance) %>%
	  select(-correlation, -distance)
 
    # using Fisher r-to-z transformation to average correlation scores
	corIDW$z <- 0.5 * (log(corIDW$filt + 1) - log(1 - corIDW$filt))	
	# average z-scores by breaks
	zAv <- ddply(corIDW, .(variable, breaks),
				summarise,
				zav = mean(z, na.rm = T))
	# back-transform the averaged z-scores
	invDistFilter <- zAv %>%
	        mutate(filter = (exp(2 * zav)-1)/(exp(2 * zav)+1)) %>%
			select(-zav)


    # add in pol columns
   # invDistFilter$pol <- str_sub(invDistFilter$site, start=4L)
	#invDistFilter$site <- str_sub(invDistFilter$site, end = 2L)


  # if alarm is in play then this section is added 
 # if(!is.na(correlation.threshold)){
 # invDistFilter <- invDistFilter %>%
#  mutate(filter = ifelse(filter < correlation.threshold, 1, 0), 
 #        test = 'idw correlation') %>%
#  filter(filter == 1)
 # } else {
#  invDistFilter$test = NA
    
 # }
  
  # return invDistFilter as the product of the function
#  invDistFilter

}

