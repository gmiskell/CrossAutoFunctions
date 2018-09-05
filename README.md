# CrossAutoFunctions
## Analysing incoming data from a network of air quality instruments.

![](https://github.com/gmiskell/CrossAutoFunctions/blob/master/auto_cross_image.png)

The package sets the data up as either: <br>
- **Cross** Where results are compared to the latest results from across the network. <br>
- **Auto** Where results are compared to the recent history from the examined location.

A collection of data-driven tests (&gamma;) are provided for checking the data, which could be weighted (&omega;) and summed to give a single value reporting on reliability. 

### Set up.

	library(devtools)
	# load from GitHub
	install_github('gmiskell/CrossAutoFunctions')
	library(CrossAutoFunctions)

Some example data, `example.RData` has been included for demonstration purposes. 

	load('/example.RData')
	example
	# A tibble: 22,362 x 5
	   date                site        lat   lon   pol
	   <dttm>              <chr>     <dbl> <dbl> <dbl>
	 1 2017-09-10 00:00:00 AQY-AA011 -43.5  173.  6.20
	 2 2017-09-11 00:00:00 AQY-AA011 -43.5  173.  3.19
	 3 2017-09-12 00:00:00 AQY-AA011 -43.5  173.  5.34
	 4 2017-09-13 00:00:00 AQY-AA011 -43.5  173.  5.56
	 5 2017-09-14 00:00:00 AQY-AA011 -43.5  173.  1.53
	 6 2017-09-15 00:00:00 AQY-AA011 -43.5  173.  3.72
	 7 2017-09-16 00:00:00 AQY-AA011 -43.5  173.  3.09
	 8 2017-09-17 00:00:00 AQY-AA011 -43.5  173.  4.28
	 9 2017-09-18 00:00:00 AQY-AA011 -43.5  173.  5.31
	10 2017-09-19 00:00:00 AQY-AA011 -43.5  173.  3.98
	# ... with 22,352 more rows

Some tests use a 'auto' and 'cross' comparison methodology (auto is where the data is compared to it's previous history and cross is where the data is compared to the network for the same time). Others have used rank-data statistics, a non-parametric method that removes the specific value of an observation and retains ordinal value in relation to the rest of the network. Following each test, they are summed by their binary output (where exceeding a threshold gives a score of 1, and 0 otherwise) and divided by the number of tests operating. This gives a proportion of how many tests are sending warnings/alarms at any time.

All tests are constructed as running functions stepping every hour, with some tests using a fixed time window size and others using daily data (sample size = 1:24). This is to try and exploit different patterns in the data and to reduce nonindependence of tests results with each other. The results from one test run provides the warning signal. This is then compared to the previous warning signals over a second fixed time window, with the sum of the binary results providing the alarm signal. This provides an idea on persistence in results, where repeated warnings can indicate systematic error. Time window sizes were determined to obtain a sample that was large enough to estimate the distribution and minimise the effect of missing values, yet small enough that response could be taken within a reasonable time-frame. Currently, the two time-windows are 72 (hours) for the warning, or test size, and 120 (hours) for the alarm, or persistence signal. 

The `autoFUN` and `crossFUN` convert the data to formats suitable for subsequent analysis.

	auto_df <- autoFUN(example, date = 'date', obs = 'pol', ell = 5, group = 'site')
	cross_df <- crossFUN(example, date = 'date', obs = 'pol', group = 'site')

### The different tests

#### The Kolmogorov-Smirnov test
This test works by using a two-sample Kolmogorov-Smirnov (KS) test, where two sets of data are compared (the observation and the proxy - here the 'network' median^1) to assess if their underlying distributions are different. The strength of using such a test for drift detection within a continuous data stream is it is robust to missing or slightly lagged data and sensitive to repetitive values, which can be assumed to become more unlikely over time due to expected natural variation in environmental monitoring. The KS test hypothesises the two samples will have similar distributions by comparing the maximum absolute distance between the y-axis values for similar x-axis values on the cumulative probability curve. KS statistics can be translated into p-values (p<sub>KS</sub>) based on sample size and probability thresholds, with the null hypothesis, H<sub>0</sub>, being the two samples could come from the same distribution^2. NOTE: it is important to remember that this test is based on a probability distribution, and so statistical significance does not necessarily translate into practical differences.

^1: 'network' here can refer to either the data across the network for the same time ('cross') or the data from previous days at the same location ('auto').

^2: here a difference between the proxy and the instrument may not always represent an error, as the proxy may not be suitable for representing concentrations at this place.

	# find the network median first
	auto_df <- auto_df %>% 
		group_by(date, site) %>% 
		do(networkmedianFUN(., obs = 'comparison.value', group = 'comparison')) %>% 
		ungroup() %>% 
		dplyr::filter(group == 'l 0')
	cross_df <- cross_df %>% 
		group_by(date) %>% 
		do(networkmedianFUN(., obs = 'comparison.value', group = 'comparison')) %>% 
		ungroup() %>% 
		rename(site = group)
	# run the ksFUN
	auto_ks <- as.data.frame(auto_df) %>% 
		group_by(site) %>% 
		do(ksFUN(., obs = 'obs', date = 'date', proxy = 'network.proxy', window.length = 72, theta = 0.05, tau = 120)) %>% 
		mutate(set.up = 'auto')
	cross_ks <- as.data.frame(cross_df) %>% 
	group_by(site) %>% 
	do(ksFUN(., obs = 'obs', date = 'date', proxy = 'network.proxy', window.length = 72, theta = 0.05, tau = 120)) %>% 
	mutate(set.up = 'cross')
	
#### The ranked variance test
The ranked variance (RV) test works by using data from all instruments, with the daily variance the variable of interest. Daily variance is a useful measure for error identification as diurnal cycles are present in most pollutants and can have similar magnitudes across an area for the same time. Ranking site variance based on their relative position across the network for each day was completed and compared on their cumulative probability curve (0 - 1 score). Those locations that either remained low (so had ongoing low diurnal variation) or those that changed in their ranking over time (so had a change in diurnal patterns that were not observed at other locations) were identified by using a threshold on the cumulative probability distribution.

	# run the rvFUN
	auto_rv <- rvFUN(auto_df, obs = 'pol', date = 'date', group = 'site', theta = 0.3, tau = 5) %>% 
		mutate(set.up = 'auto')
	cross_rv <- rvFUN(cross_df, obs = 'pol', date = 'date', group = 'site', theta = 0.3, tau = 5) %>% 
		mutate(set.up = 'cross')

#### The percentile test
The percentile (P) test works by using data from all instruments, with the daily Median Absolute Deviation (MAD) score used. This is similar to a z-score, or effect size, and so illustrates those locations that have practically significant differences. The daily minima and maxima were calculated for each location and their MAD score calculated for each, relative to the rest of the network. This would illustrate where locations had differences in their lower and higher values, in order to capture where drifting baselines, or high concentrations were occurring. Here, differences greater than 2&sigma; were used as a threshold, with two outputs (the minima and maxima score).

	# run the percentFUN
	auto_p <- percentFUN(auto_df, obs = 'pol', date = 'date', group = 'site', theta = 2, tau = 5) %>% 
		mutate(set.up = 'auto')
	cross_p <- percentFUN(cross_df, obs = 'pol', date = 'date', group = 'site', theta = 2, tau = 5) %>% 
		mutate(set.up = 'cross')
		
#### The decomposed regression test
The decomposed regression (DR) test works by removing an estimate on the long- and short-term trends in the data and looking at how the residuals, or error terms, change between the instrument and the proxy (here the 'network' median^1). Following the calculation of these error terms, where we can expect to see occasional random error due to local variability, we remove one error term from the other (here instrument &epsilon; - proxy &epsilon;) to find where they are different and especially where they remain different^2. The long-term trend used here is the site's weekly average and the short-term trend used here is the daily autocorrelation function (n = 24 when using hourly averages and so 24 coefficients in the ACF).

^1: 'network' here can refer to either the data across the network for the same time ('cross') or the data from previous days at the same location ('auto'). 

^2: here a difference between the proxy and the instrument may not always represent an error, as the proxy may not be suitable for representing concentrations at this place.

	# run the dfFUN
	# find the network median, as in the ksFUN example
	auto_dr <- as.data.frame(auto_df) %>% 
		group_by(site) %>% 
		do(drFUN(., obs = 'obs', date = 'date', proxy = 'network.proxy', long.term = 168, short.term = 24, theta = 40, tau = 120)) %>% 
		mutate(set.up = 'auto')
	cross_dr <- as.data.frame(cross_df) %>% 
		group_by(site) %>% 
		do(drFUN(., obs = 'obs', date = 'date', proxy = 'network.proxy', long.term = 168, short.term = 24, theta = 40, tau = 120)) %>% 
		mutate(set.up = 'cross')
		
Results are set up so that they can be combined into one data frame with a mean of the number of alarms per *t*.

	all_data <- bind_rows(auto_ks, cross_ks, auto_rv, cross_rv, auto_p, cross_p, auto_dr, cross_dr) %>% 
		group_by(date, site) %>% 
		mutate(warning_sum = mean(warning, na.rm = T), alarm_sum = mean(alarm, na.rm = T) %>% 
		ungroup()