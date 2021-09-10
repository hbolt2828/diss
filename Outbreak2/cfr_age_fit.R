####################################################################################
##                                                                                ##
##  2) Model Fitting and Parameter Inference                                      ##
##                                                                                ##
##    The age-disaggregated incidence data for NorthWest and outside is then    	##
##    integrated with age-disaggregated deaths from NorthWest and outside over 	  ##
##    the same time period. After a age specific CFR model is fitted to find how  ##
##    the CFR varies. This must take into account under reporting of cases and    ##
##    age-specific patterns of disease severity                                   ##
##                                                                                ##
####################################################################################

# Installing packages and prior and likelihood 
set.seed(101092)
source("~/desktop/mcmc_with_my_prior.R")

library(devtools)
devtools::install_github("mrc-ide/drjacoby", ref = "version1.0")  # uncomment and run this line the first time through
library(drjacoby)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Loading in Age and Location Disaggregated Case Data
age_disaggregated_case_onset <- readRDS("~/desktop/age_disaggregated_onset_incidence_data.rds")
age_cases <- age_disaggregated_case_onset %>%
  group_by(age_groups) %>%
  summarise(counts = sum(cases))
data <- age_disaggregated_case_onset %>%
  arrange(location, age_groups, date)

# Setting Static Variables
age_groups <- c("0_9", "10_19", "20_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80+")
n_age_bands <- length(unique(data$age_groups))
max_date <- as.numeric(max(data$date) - min(data$date)) # most recent date in the date, relative to date start
observed_deaths <- c(0, 5, 19, 55, 174, 471, 1243, 3269, 7714) # (by age of UK, summing up England, Welsh, NI, Scotland data from Gov. UK
prop_deaths_NorthWest <- 3632 / 12950  # Using sum of values above and North West dates over this period from Gov.UK. #15492
# apply this proportion to the number of deaths detailed in the UK.
deaths_NorthWest <- round(prop_deaths_NorthWest * sum(observed_deaths))
deaths_outside <- (1 - prop_deaths_NorthWest) * sum(observed_deaths)
death_observation_censoring <- as.numeric(as.Date("2020-01-29") - min(data$date)) # Assume no deaths detecte before 29th of January 2020 in the UK 
deaths_x <- round(prop_deaths_NorthWest * sum(observed_deaths))
deaths_n <- sum(observed_deaths)
# Using results of deaths and cases from Gov.UK data:
Gov_observed_deaths <- 12950 
Gov_observed_cases <- 1015927

# Input data in our data frame
x <- c(0, n_age_bands,
       death_observation_censoring,
       min(data[data$location == "Outside", "nici"]),
       deaths_x, deaths_n,
       Gov_observed_deaths, Gov_observed_cases,
       observed_deaths,
       data$cases, as.numeric(as.factor(data$age_groups)),
       data$date - min(data$date), as.numeric(data$location == "NorthWest"),
       data$nici)

# Input data for output = prediction
x_output <- x
x_output[1] <- 1
saveRDS(x_output, "~/desktop/predicted_x.rds")

#   Setting our parameters min and max values for our MCMC to find the best value 
#   z scales NorthWest relative to outside
#   r growth rate (fixed)
#   D detection window (fixed)
df_paramsRR <- data.frame(name = c("m_od", "s_od", "maxday", "z", "r", "D", "cfr_80+", "RR_0_9", "RR_10_19", "RR_20_29", "RR_30_39", "RR_40_49", "RR_50_59", "RR_60_69", "RR_70_79"),
                          min = c(10, 0, max_date, 0, 0, 0, rep(0, 9)),
                          max = c(Inf, Inf, max_date, 1, 0.1, 14, 1, rep(Inf, 8)),
                          init = c(15, 0.35, max_date, 0.005, 0.05, 10, 0.1, rep(1, 8)))
params <- df_paramsRR$init
lL_plane(df_paramsRR$init, x)

# Using MCMC, stating how many samples need to be taken and how many we will use  
burnin <- 2000
samples <- 1000
mcmc_output <- run_mcmc(data = x,
                        df_params = df_paramsRR,
                        loglike = lL_plane,
                        logprior = lP_plane,
                        burnin = burnin,
                        samples = samples,
                        chains = 1)
# (Run with chains > 1 for convergence checks)
saveRDS(mcmc_output, "~/desktop/complete_output.rds")
saveRDS(mcmc_output$output, "~/desktop/MCMC_fitting_output.rds") # save the MCMC output
saveRDS(list(inputs = x, 
             parameters = df_paramsRR,
             samples = samples,
             burnin = burnin), "~/desktop/MCMC_inputs_and_parameters.rds") # save the data list used to run the MCMC

