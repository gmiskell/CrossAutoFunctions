#' A daily ranked variance function
#'
#' This function gives the daily ranked variance for an observation, relative to a pool of other observations.
#' @param reflective This decides whether to use a selected time period of data (reflective), or to use data from the latest day (live). Defaults to TRUE.
#' @param plot This decides whether to plot the data. Defaults to TRUE.
#` @export
#' @examples
#' ARFUN()

# AR, MR script using MA work

# AR function

# x = data file where the data of interest are in a column - no differentiation between sites
#     and so this needs to be done outside of the function
# dat = column of interest in data file x
# M = length of training data
# N = length of autocorrelation function + 1 (for the prediction column)


ARFUN <- function(x, dat = 'conc_mean', M = 168, N = 25){
 
  # load required libraries
  library(stats); library(raster); library(dplyr)
  
  # training data - select first [M] length of data
  
  z.train <- x %>%
    head(n = M)
  
  # derive the average
  
  z.train.av <- mean(z.train[, dat], na.rm = T)
  
  # remove the average from the observations & force NA to 0 for matrix calculations
  
  x.obs <- z.train[, dat] - z.train.av
  x.obs[is.na(x.obs)] <- 0
  
  # short-term autocorrelation matrix to [N] length
  
  x.obs.lag <- embed(x.obs, N)
  
  # X here is the data to be predicted, H here is the data used for predictions
  x.obs.X <- x.obs.lag[, 1]
  x.obs.H <- x.obs.lag[, N:2]
  
  # matrix coefficients (N-1 length)
  b <- solve(qr(x.obs.H), x.obs.X)
   
  
  # application of above to the following data
  
  # derive M length averages - these will be rolling at each unit
  # the first M data will use the set average
  
  df <- x %>%
    mutate(movingav = movingFun(x[, dat], M, fun = mean, type = 'to', na.rm = T)) %>%
    mutate(movingav = replace(movingav, 1:M, z.train.av)) %>%
    mutate(X = x[, dat] - movingav)
  
  # set up the matrix
  x.obs.lag <- embed(df$X, N)
  x.obs.X <- x.obs.lag[, 1]
  x.obs.H <- x.obs.lag[, N:2]
  
  # solve the matrix using the derived coefficients from the model
  b.X <- sweep(x.obs.H, MARGIN = 2, b, '*')
  # this sums all coefficient predictions at each unit to give an estimate of column N+1
  x.Hb <- rowSums(b.X, na.rm = T)
  # determine the difference between observed and predicted
  iX <- x.obs.X - x.Hb
  # insert some leading N values
  rep.NA <- rep(NA, N-1)
  iX <- c(rep.NA, iX)
  # join results to original data file
  iX.data <- cbind(df, iX) 
  
  # print trace
  plot(iX.data$iX, type = 'l')
    
}


