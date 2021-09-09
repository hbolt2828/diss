# Load packages
library(tidyverse)
library(forecast)
library(smoother)
#Putting in pre sets for gamma value for Bettencourt and Ribeiros formula
# Setting max r value, realistically will rarely be above 2, but 12 is used for safety

# r_t_range is a vector of possible values for R_t 
R_T_MAX = 12
r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)

# Gamma is 1/serial interval as un Du 's paper
GAMMA = 1/4

# Inserting English data from Gov.UK for regional rates
 covid_infections <- read.csv(file = 'data/regional_r_rates.csv')
 covid_infections <- readr::read_csv('data/regional_r_rates.csv')
 head(covid_infections)
 

 # Compute new infections and smooth them as the daily rates may vary for under ascertainment from 
 smooth_new_infections <- function(infections){
   infections %>%
     arrange(date) %>%
     #find the number of new daily infections by subtracting from the previous day
     mutate(new_infections = c(infections[1], diff(infections))) %>%
     mutate(new_infections_smooth = round(
       smoother::smth(new_infections, window = 7, tails = TRUE)
     )) %>%
     select(region, date, new_infections, new_infections_smooth)
 }
 
 #selecting North West as an example region
 region_selected <- "North West"
 #seperate the data to only give that of the North West and then smooth
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   head()
 
 
 # Plotting number of new infections in the North West region showing the daily trend
 # the dotted line is the true number of infections recorder per day, whilst the connected line is the smoothed data
 plot_new_infections <- function(infections){
   infections %>%
     ggplot(aes(x = date, y = new_infections, group = 1)) +
     geom_line(linetype = 'dotted', color = 'gray40') +
     geom_line(aes(y = new_infections_smooth), color = "#14243e") +
     labs(
       title = "New infections per day",
       subtitle = unique(infections$region),
       x = NULL, y = NULL
     )
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   plot_new_infections()
 
 
 #Finding the likelihood for each time period, using the Poisson distribution
 compute_likelihood <- function(infections){
   likelihood <- infections %>%
     filter(new_infections_smooth > 0) %>%
     mutate(
       r_t = list(r_t_range),
      #Bettencourt and Ribeiro's formula for lambda
       lambda = map(lag(new_infections_smooth, 1), ~ .x * exp(GAMMA * (r_t_range - 1))),
       likelihood_r_t = map2(new_infections_smooth, lambda, dpois, log = TRUE)
     ) %>%
     slice(-1) %>%
     select(-lambda) %>%
     unnest(c(likelihood_r_t, r_t))
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   compute_likelihood() %>%
   head()
 
 
 #Computing the posterior 
 compute_posterior <- function(likelihood){
   likelihood %>%
     arrange(date) %>%
     group_by(r_t) %>%
     mutate(posterior = exp(
       zoo::rollapplyr(likelihood_r_t, 7, sum, partial = TRUE)
     )) %>%
     group_by(date) %>%
     mutate(posterior = posterior / sum(posterior, na.rm = TRUE)) %>%
     # HACK: NaNs in the posterior create issues later on. So we remove them.
     mutate(posterior = ifelse(is.nan(posterior), 0, posterior)) %>%
     ungroup() %>%
     select(-likelihood_r_t)
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   compute_likelihood() %>%
   compute_posterior() %>%
   head()
 
 
 #Plotting the posterior 
 plot_posteriors <- function(posteriors){
   posteriors %>%
     ggplot(aes(x = r_t, y = posterior, group = date)) +
     geom_line(alpha = 0.2) +
     labs(
       title = expression(paste("Daily Posterior of R"[t], " by day")),
       subtitle = unique(posteriors$region),
       x = '',
       y = ''
     ) +
     coord_cartesian(xlim = c(0.4, 4)) +
     theme(legend.position = 'none')
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   compute_likelihood() %>%
   compute_posterior() %>%
   plot_posteriors()
 
 # We can see the result from the posterior that the mean is around 1. 
 # Estimating the value for R_t,  and using 10000 samples from the posterior to find the 95% HDPI, large
 estimate_rt <- function(posteriors){
   posteriors %>%
     group_by(region, date) %>%
     summarize(
       r_t_simulated = list(sample(r_t_range, 10000, replace = TRUE, prob = posterior)),
       r_t_most_likely = r_t_range[which.max(posterior)]
     ) %>%
     mutate(
       r_t_lo = map_dbl(r_t_simulated, ~ hdi(.x)[1]),
       r_t_hi = map_dbl(r_t_simulated, ~ hdi(.x)[2])
     ) %>%
     select(-r_t_simulated)
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   compute_likelihood() %>%
   compute_posterior() %>%
   estimate_rt() %>%
   head()
 
 
 #Ploting our estimates
 plot_estimates <- function(estimates){
   estimates %>%
     ggplot(aes(x = date, y = r_t_most_likely, group = 1)) +
     geom_point(color = "darkorange", alpha = 0.8, size = 0.5) +
     geom_line(color = "#14243e") +
     geom_hline(yintercept = 1, linetype = 'dashed') +
     geom_ribbon(
       aes(ymin = r_t_lo, ymax = r_t_hi),
       fill = 'darkred',
       alpha = 0.2
     ) +
     labs(
       title = expression('Real time R'[t]), x = '', y = '',
       subtitle = unique(estimates$region)
     ) +
     coord_cartesian(ylim = c(0, 4))
 }
 
 covid_infections %>%
   filter(region == region_selected) %>%
   smooth_new_infections() %>%
   compute_likelihood() %>%
   compute_posterior() %>%
   estimate_rt() %>%
   plot_estimates()
 

 
 
 # This function can take a couple of minutes to run
 #   as it loops across all regions
 estimates_all <- covid_infections %>%
   filter(date >= "2020-03-21") %>%
   group_by(region) %>%
   # Ignore regions that have not reached 100 infections
   filter(max(infections) > 100 ) %>%
   group_split() %>%
   map_df(~ {
     .x %>%
       smooth_new_infections() %>%
       compute_likelihood() %>%
       compute_posterior() %>%
       estimate_rt()
   }) %>%
   ungroup()
 
 estimates_all %>%
   head()
 
 
 # Increase plot height and width
 options(repr.plot.height = 40, repr.plot.width = 20)
 estimates_all %>%
   plot_estimates() +
   facet_wrap(~ region, ncol = 3) +
   labs(subtitle = "")
 
 # Reset plot dimensions
 options(repr.plot.height = 12, repr.plot.width = 8)
 
library(geofacet)
 #⚠️ CAUTION: This code will error on colab
options(repr.plot.height = 40, repr.plot.width = 20)
estimates_all %>%
mutate(region = region.abb[match(region, region.name)]) %>%
plot_estimates() +
geofacet::facet_geo(~ region, ncol = 4) +
labs(subtitle = "") +
theme(strip.text = element_text(hjust = 0))
options(repr.plot.height = 12, repr.plot.width = 8)
 
 options(repr.plot.width = 20, repr.plot.height = 8)
 estimates_all %>%
   group_by(region) %>%
   filter(date == max(date)) %>%
   ungroup() %>%
   mutate(region = forcats::fct_reorder(region, r_t_most_likely)) %>%
   ggplot(aes(x = region, y = r_t_most_likely))  +
   geom_hline(yintercept = 1, linetype = 'dotted') +
   geom_errorbar(aes(ymin = r_t_lo, ymax = r_t_hi), width = 0.2) +
   scale_fill_manual(values = c(None = 'darkred', Partial = 'gray50', Full = 'gray70')) +
   labs(
     title = expression(paste("Most Recent R"[t], " by region")),
     x = '', y = ''
   ) +
   theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5))
 options(repr.plot.width = 12, repr.plot.height = 5)
 
 
 
 

 ########################Future Predictions############################
 
 
 #east mids
 
 east_mids_r <- estimates_all[ which(estimates_all$region=='East Midlands'), ]
 east_mids_r1 <- data.frame(y = 1:504, east_mids_r)
 east_mids_r2 <-subset(east_mids_r1, y <= 504)
 east_mids_r3 <-subset(east_mids_r1, y > 504)
 
 
 
 
#Plot time series for all records, infections and deaths
 
#Auto arima to forecast the future of all 9 regions
 
 ts_east_mids_r2 <-ts (east_mids_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_east_mids_r2)
 #Registered infections forcasting
 autoarima1 <- auto.arima(ts_east_mids_r2)
 forecast1 <- forecast(autoarima1, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)

 plot(forecast1,ylab='R rate over time in the East Midlands',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 #east of england
 
 east_england_r <- estimates_all[ which(estimates_all$region=='East of England'), ]
 east_england_r1 <- data.frame(y = 1:505, east_england_r)
 east_england_r2 <-subset(east_england_r1, y <= 505)
 east_england_r3 <-subset(east_england_r1, y > 505)
 
 #Plot time series for all records, infections and deaths#
 ts_east_england_r2 <-ts (east_england_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_east_england_r2)
 #Registered infections forcasting
 autoarima2 <- auto.arima(ts_east_england_r2)
 forecast2 <- forecast(autoarima2, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast2,ylab='R rate over time in the East of England',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 
 
 
 #Ldn
 
 ldn_r <- estimates_all[ which(estimates_all$region=='London'), ]
 ldn_r1 <- data.frame(y = 1:504, ldn_r)
 ldn_r2 <-subset(ldn_r1, y <= 504)
 ldn_r3 <-subset(ldn_r1, y > 504)
 
 #Plot time series for all records, infections and deaths#
 ts_ldn_r2 <-ts (ldn_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_ldn_r2)
 #Registered infections forcasting
 autoarima3 <- auto.arima(ts_ldn_r2)
 forecast3 <- forecast(autoarima3, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)

 plot(forecast3,ylab='R rate over time in the London',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 
 
 
 #North East
 
 NE_r <- estimates_all[ which(estimates_all$region=='North East'), ]
 NE_r1 <- data.frame(y = 1:504, NE_r)
 NE_r2 <-subset(NE_r1, y <= 504)
 NE_r3 <-subset(NE_r1, y > 504)
 
 #Plot time series for all records, infections and deaths#
 ts_NE_r2 <-ts (NE_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_NE_r2)
 #Registered infections forcasting
 autoarima4 <- auto.arima(ts_NE_r2)
 forecast4 <- forecast(autoarima4, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast4,ylab='R rate over time in the North East',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 #North West
 
 NW_r <- estimates_all[ which(estimates_all$region=='North West'), ]
 NW_r1 <- data.frame(y = 1:504, NW_r)
 NW_r2 <-subset(NW_r1, y <= 504)
 NW_r3 <-subset(NW_r1, y > 504)
 
 #Plot time series for all records, infections and deaths#
 ts_NW_r2 <-ts (NW_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_NW_r2)
 #Registered infections forcasting
 autoarima5 <- auto.arima(ts_NW_r2)
 forecast5 <- forecast(autoarima5, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast5,ylab='R rate over time in the North West',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 #South East
 
 SE_r <- estimates_all[ which(estimates_all$region=='South East'), ]
 SE_r1 <- data.frame(y = 1:505, SE_r)
 SE_r2 <-subset(SE_r1, y <= 505)
 SE_r3 <-subset(SE_r1, y > 505)
 
 #Plot time series for all records, infections and deaths#
 ts_SE_r2 <-ts (SE_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_SE_r2)
 #Registered infections forcasting
 autoarima6 <- auto.arima(ts_SE_r2)
 forecast6 <- forecast(autoarima6, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast6,ylab='R rate over time in the South East',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 
 #South West
 
 SW_r <- estimates_all[ which(estimates_all$region=='South West'), ]
 SW_r1 <- data.frame(y = 1:505, SW_r)
 SW_r2 <-subset(SW_r1, y <= 505)
 SW_r3 <-subset(SW_r1, y > 505)
 
 #Plot time series for all records, infections and deaths#
 ts_SW_r2 <-ts (SW_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_SW_r2)
 #Registered infections forcasting
 autoarima7 <- auto.arima(ts_SW_r2)
 forecast7 <- forecast(autoarima7, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast7,ylab='R rate over time in the South West',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 # West Midlands
 
 WM_r <- estimates_all[ which(estimates_all$region=='West Midlands'), ]
 WM_r1 <- data.frame(y = 1:505, WM_r)
 WM_r2 <-subset(WM_r1, y <= 505)
 WM_r3 <-subset(WM_r1, y > 505)
 
 #Plot time series for all records, infections and deaths#
 ts_WM_r2 <-ts (WM_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_WM_r2)
 #Registered infections forcasting
 autoarima8 <- auto.arima(ts_WM_r2)
 forecast8 <- forecast(autoarima8, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast8,ylab='R rate over time in the West Midlands',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 # Yorkshire and the Humber
 
 y_r <- estimates_all[ which(estimates_all$region=='Yorkshire and The Humber'), ]
 y_r1 <- data.frame(y = 1:503, y_r)
 y_r2 <-subset(y_r1, y <= 503)
 y_r3 <-subset(y_r1, y > 503)
 
 #Plot time series for all records, infections and deaths#
 ts_y_r2 <-ts (y_r2$r_t_most_likely, frequency = 1,start = c(22/03/2020,1))
 plot(ts_y_r2)
 #Registered infections forcasting
 autoarima9 <- auto.arima(ts_y_r2)
 forecast9 <- forecast(autoarima9, h=30) # pretty sure h is the number of future points forecasted
 #calculate MAPE (mean absolute percentage error)
 
 plot(forecast9,ylab='R rate over time in the Yorkshire and the Humber',xlab='Days from 2020-03-22')
 abline(h = 1, col='purple')
 
 par(mfrow=c(3,3))
 plot(forecast1,ylab='East Midlands',xlab='Days from 2020-03-22')
 plot(forecast2,ylab='East of England',xlab='Days from 2020-03-22')
 plot(forecast3,ylab='London',xlab='Days from 2020-03-22')
 plot(forecast4,ylab='North East',xlab='Days from 2020-03-22')
 plot(forecast5,ylab='North West',xlab='Days from 2020-03-22')
 plot(forecast6,ylab='South East',xlab='Days from 2020-03-22')
 plot(forecast7,ylab='South West',xlab='Days from 2020-03-22')
 plot(forecast8,ylab='West Midlands',xlab='Days from 2020-03-22')
 plot(forecast9,ylab='Yorkshire and the Humber',xlab='Days from 2020-03-22')
 
 
 
# predicting future using sampling
 
 #East Mids
 drawsDF1 = east_mids_r
 randomDrawDF1 = drawsDF1 %>% 
   sample_n(size = 50)
 approx1 <- data.frame(y = 504:553, randomDrawDF1)
 plot(ts_east_mids_r2)
 lines(approx1$y, approx1$r_t_most_likely, type = "l", col = "blue")
 
 
 #East of Eng
 drawsDF2 = east_england_r
 randomDrawDF2 = drawsDF2 %>% 
   sample_n(size = 50)
 approx2 <- data.frame(y = 504:553, randomDrawDF2)
 plot(ts_east_england_r2)
 lines(approx2$y, approx2$r_t_most_likely, type = "l", col = "blue")
 
 
 #London
 drawsDF3 = ldn_r
 randomDrawDF3 = drawsDF3 %>% 
   sample_n(size = 50)
 approx3 <- data.frame(y = 504:553, randomDrawDF3)
 plot(ts_ldn_r2)
 lines(approx3$y, approx3$r_t_most_likely, type = "l", col = "blue")
 
 
 #North East
 drawsDF4 = NE_r
 randomDrawDF4 = drawsDF4 %>% 
   sample_n(size = 50)
 approx4 <- data.frame(y = 504:553, randomDrawDF4)
 plot(ts_NE_r2)
 lines(approx4$y, approx4$r_t_most_likely, type = "l", col = "blue")
 
 #North West
 drawsDF5 = NW_r
 randomDrawDF5 = drawsDF5 %>% 
   sample_n(size = 50)
 approx5 <- data.frame(y = 504:553, randomDrawDF5)
 plot(ts_NW_r2)
 lines(approx5$y, approx5$r_t_most_likely, type = "l", col = "blue")
 
 #South East
 drawsDF6 = SE_r
 randomDrawDF6 = drawsDF6 %>% 
   sample_n(size = 50)
 approx6 <- data.frame(y = 504:553, randomDrawDF6)
 plot(ts_SE_r2)
 lines(approx6$y, approx6$r_t_most_likely, type = "l", col = "blue")

 
  #South West 
 drawsDF7 = SW_r
 randomDrawDF7 = drawsDF7 %>% 
   sample_n(size = 50)
 approx7 <- data.frame(y = 504:553, randomDrawDF7)
 plot(ts_SW_r2)
 lines(approx7$y, approx7$r_t_most_likely, type = "l", col = "blue")
 
 
 #West Midlands
 drawsDF8 = WM_r
 randomDrawDF8 = drawsDF8 %>% 
   sample_n(size = 50)
 approx8 <- data.frame(y = 504:553, randomDrawDF8)
 plot(ts_WM_r2)
 lines(approx8$y, approx8$r_t_most_likely, type = "l", col = "blue")
 
 
 
 #Yorkshire and the humber
 
 drawsDF9 = y_r
 randomDrawDF9 = drawsDF9 %>% 
   sample_n(size = 50)
 approx9 <- data.frame(y = 504:553, randomDrawDF9)
 
 plot(ts_y_r2)
 lines(approx9$y, approx9$r_t_most_likely, type = "l", col = "blue")
 
 
 #Plotting a 3x3 of the time series with the new future approximations
 
 par(mfrow=c(3,3))
 
 plot(ts_east_mids_r2, ylab = "East Midlands R rate", xlab = "days since 2020-03-22")
 lines(approx1$y, approx1$r_t_most_likely, type = "l", col = "blue")
 
 plot(ts_east_england_r2, ylab = "East of England R rate", xlab = "days since 2020-03-22")
 lines(approx2$y, approx2$r_t_most_likely, type = "l", col = "blue")

  plot(ts_ldn_r2, ylab = "London R rate", xlab = "days since 2020-03-22")
 lines(approx3$y, approx3$r_t_most_likely, type = "l", col = "blue")
 
 plot(ts_NE_r2, ylab = "North East R rate", xlab = "days since 2020-03-22")
 lines(approx4$y, approx4$r_t_most_likely, type = "l", col = "blue")

 plot(ts_NW_r2, ylab = "North West R rate", xlab = "days since 2020-03-22")
 lines(approx5$y, approx5$r_t_most_likely, type = "l", col = "blue")
 
 plot(ts_SE_r2, ylab = "South East R rate", xlab = "days since 2020-03-22")
 lines(approx6$y, approx6$r_t_most_likely, type = "l", col = "blue")
 
 plot(ts_SW_r2, ylab = "South West R rate", xlab = "days since 2020-03-22")
 lines(approx7$y, approx7$r_t_most_likely, type = "l", col = "blue")

 plot(ts_WM_r2, ylab = "West Midlands R rate", xlab = "days since 2020-03-22")
 lines(approx8$y, approx8$r_t_most_likely, type = "l", col = "blue")
 
 plot(ts_y_r2, ylab = "Yorkshire and the Humber R rate", xlab = "days since 2020-03-22")
 lines(approx9$y, approx9$r_t_most_likely, type = "l", col = "blue")
 
