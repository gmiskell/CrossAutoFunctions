#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' inverseDistanceFilter()
# inverse distance correlation filter script

# description:  this filter will compare the correlations of a site against the other 
#               sites. a weighting will be applied, giving higher value to those sites 
#               closer. the final product will be a sum of the correlations multiplied 
#               by the inverse distance

inverseDistanceFilter <- function(df, meta = meta.data, correlation.threshold = NA, 
                                  window.length = 1440){
  
  # load libraries
  library(openair);library(plyr);library(stringr);
  library(zoo);library(reshape2);library(dplyr)

  # rearrange the data to perform correlations per column
  idw.filter <- df %>%
    select(date, site, pol, conc_mean) %>%
    dcast(... ~ site + pol, median)
  

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
    select(-value) %>%
    mutate(variable = str_c(site, '.', variable))


  # Spearman rank correlations
  date <- idw.filter$date
  # allow for a lag in the data for rolling functions
  date <- date[-c(1:(window.length-1))]
  # convert the data into a zoo, a class in R used for rolling functions
  z <- zoo(idw.filter[,-1])
  # gather names of the correlations
  c <- cor(z, use='pairwise.complete.obs')
  l <- as.logical(c)
  n <- paste(rownames(c)[row(c)[l]], rownames(c)[col(c)[l]])
  # perform correlation on concentrations using data where it is recorded at both locations
  roll.corFUN <- rollapply(z, width = window.length, 
                           FUN = function(Z) { 
                             return(cor(Z, use = 'pairwise.complete.obs', 
                                        method = 'spearman')[l])}, 
                           by.column = F, align = 'right')
  # add names and bind to date column
  colnames(roll.corFUN) <- n
  d <- data.frame(date, roll.corFUN)


  # join distance and correlations and create filter
  corIDW <- d %>%
    melt(id.vars = 'date') %>%
    rename(correlation = value) %>%
    join(dist.melt, by = 'variable') %>%
    # remove correlations performed on the same site and NA correlations
    filter(distance > 0 & !is.na(correlation)) %>%
    # add in the inverse distance
    mutate(filt = correlation * 1 / distance)
 
  # sum the data into one figure
  invDistFilter <- ddply(corIDW, .(site, date),
                       summarise,
                       filter = sum(filt, na.rm = T))

  # add in pol and site columns
  invDistFilter$pol <- str_sub(invDistFilter$site, start=4L)
  invDistFilter$site <- str_sub(invDistFilter$site, end=2L)

  # if alarm is in play then this section is added 
  if(!is.na(correlation.threshold)){
  invDistFilter <- invDistFilter %>%
  mutate(filter = ifelse(filter < correlation.threshold, 1, 0), 
         test = 'idw correlation') %>%
  filter(filter == 1)
  } else {
  invDistFilter$test = NA
    
  }
  
  # return invDistFilter as the product of the function
  invDistFilter

}

