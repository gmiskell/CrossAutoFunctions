#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' ksFUN()
# description: rolling KS-test filter script
#              this filter compares location data to the network distibution
#              using a non-parametric statistic and prxies
#              the proxy example here is for the average for that time across the network 
#              for each pollutant (removing site under observation)


# need to create some prior functions
# proxy - create a variable network average which removes the location under interest
proxyFUN <- function(x){
  
  # load required libraries
  library(dplyr)
  
  # filter data to that of interest, removingg empty rows
  z = x %>%
    select(date, conc_mean, site, pol)
  z = na.omit(z)
  
  # make comparison data
  proxy <- unlist(sapply(seq(1, nrow(z)), function(i) {
    
    site = z[i,]
    test = z[-i,]
    
    proxy = median(test$conc_mean, na.rm = T)
    
  }))
  
  # join proxy variable to data and return
  b <- cbind(z, proxy)
  
}

# rolling KS function
rollingKStest <- function(z, var.one = 'proxy', var.two = 'conc_mean', window.1 = 72) {
  
  # load required libraries
  library(zoo); library(plyr)
    
  # interpolate data
  z[, var.one] <- na.approx(z[,var.one], na.rm = F)
  z[, var.two] <- na.approx(z[,var.two], na.rm = F)
  
  # save date column
  date <- z[, 'date']
  
  y <- data.frame(z[,var.one], z[,var.two])
  
  # convert df into zoo classes used in R for rolling functions
  y.zoo <-zoo(y)
  colnames(y.zoo) <- c('r', 's')
  
  # running ks functions
  p.value <- rollapply(y.zoo, function (x) ks.test(x[, 'r'], x[, 's'])$p.value, 
	by.column = F, fill = NA, align = 'right', width = window.1)
  statistic <- rollapply(y.zoo, function (x) ks.test(x[, 'r'], x[, 's'])$statistic, 
	by.column = F, fill = NA, align = 'right', width = window.1)
  
  p.value <- data.frame(p.value)
  statistic <- data.frame(statistic)
  
  # make results into df and tidy
  model <- cbind(p.value, statistic)
  names(model) <- c('p.value', 'statistic')
  model <- data.frame(lapply(model, function(x) round(x, 6)))
  model <- cbind(date, model)
  
  # join ks results to the data and return
  z <- join(z, model, by = 'date')
  
}


# ks function to use as final function
ks.function <- function(y, var.one = 'conc_mean', var.two = 'proxy', window.2 = 120){
  
  # load required libraries
  library(raster); library(plyr); library(stringr): library(dplyr)
  
  y <- y %>%
	select(date, site, pol, conc_mean)
	
  # generate proxy variable using function
  y <- ddply(na.omit(y), .(date, pol),
             function(x) proxyFUN(x))
  
  # generate rolling ks-test results using function by site
  ksAlarm1 <- ddply(y, .(site),
                    function(x) rollingKStest(x, var.one = var.one, var.two = var.two))
  
  # filter to those with p-values < 0.05 (significantly different to the network proxy)
  ksAlarm1$alarm1 <- ifelse(ksAlarm1$p.value < 0.05, 1, 0)
  
  #alarm generated is for systematic or local change between obs and proxy using window.2
  ksAlarm2 <- ddply(ksAlarm1, .(site, pol),
                  transform,
                  filter = movingFun(alarm1, n = window.2, fun = 'mean', type = 'to'))
  
  # filter results to those that exceed tau by 95% of recent obs.
  ksAlarm2 <- ksAlarm2 %>%
  filter(filter > 0.95)
  
  # gather variables of interest and return
  ks.filter <- ksAlarm2 %>%
  mutate(day = str_sub(date, end = -9L), test = 'ks filter') %>%
  select(day, site, pol, test, filter) %>%
  unique()
  
  ks.filter
  
}
  