  # clause on either auto or cross test to be run
  if (auto == TRUE){
  cast.dat <- x %>%
    filter(date >= as.POSIXct(date.start) - days(ell), date < (date.end)) %>%
	mutate(day = str_sub(date, end = 10)) %>%
    group_by(site) %>%
	filter(!is.na(pol)) 	
	len <- as.numeric(difftime(cast.dat$date[2],cast.dat$date[1],units = 'mins'))
	if(len == 60) len = 24
	if(len == 1) len = 1440	
	nn <- seq(1:ell) * len		
	cast.dat <- setDT(all.data.melt.hr)[, paste("l", 1:ell) := shift(value, n =  nn),by = site]
	cast.dat <- data.frame(cast.dat)
  } else {
  # rearrange data structure and create day variable
  cast.dat <- x %>%  
    filter(date >= (date.start) & date < (date.end)) %>%
    dcast(... ~ site, value.var = 'value', median) %>%
    mutate(day = str_sub(date, end = 10L))
  }
  
  
  # sort data into auto function ----

autoFN <- function(x, ell = 3, t = '2016-12-01'){
  
  x$day <- str_sub(x$date, end = 10)

  response <- x %>%
    filter(day == t) %>%
    mutate(time = str_sub(date, start = 12)) %>%
    select(time, site, resp = value)
  
  comparison <- x %>%
    filter(as.Date(day) >= as.Date(t)-ell & as.Date(day) <= as.Date(t)-1) %>%
    mutate(time = str_sub(date, start = 12)) %>%
    select(time, site, comp = value)
  
  z <- merge(response, comparison, by = 'time')
  
}

crossFN <- function(x, site2 = 'B1', t = '2016-12-01'){
  
  x$day <- str_sub(x$date, end = 10)
  
  response <- x %>%
    filter(day == t) %>%
    mutate(time = str_sub(date, start = 12)) %>%
    select(time, resp.site = site, resp = value)
  
  comparison <- x %>%
    filter(day == t & site != site2) %>%
    mutate(time = str_sub(date, start = 12)) %>%
    select(time, comp.site = site, comp = value)
  
  z <- merge(response, comparison, by = 'time')

}