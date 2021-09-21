####################################################################################
##                 Vaccinations period                                            ##
##  1) Data Collation and Preprocessing - An Overview                             ##
##                                                                                ##
##   This section loads in data on the demographic structure of the UK            ##
##   population, as well as daily case data for NorthWest (Region) and rest of    ##
##   UK from information from the UK Health Deparment. This is integrated         ##
##   with information on the Age distribution of cases from reports from the      ##
##   UK Govt to calculate the age distribution of cases for NorthWest 	          ##
##   and the rest of UK.                                                          ##
##                                                                                ##
##   Assuming an invariant age-distribution of cases over time, this is then      ##
##   used to convert the aggregate daily incidence of cases (by symptom onset)    ##
##   into age-specific daily case incidences.                                     ##
##                                                                                ##
####################################################################################

# Loading Libraries
library(tidyverse)

# Loading in Demographic and Case Data
population <- read.csv("~/desktop/UK_population_demography.csv") 
Epicentro_report_onset <- read.csv("~/desktop/Gov_UK_data.2.v.csv")

# Setting Static Variables
smoothing_centre <- as.Date("2021-06-28", format = "%Y-%m-%d")
days_to_smooth_either_side <- 3
age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

####################################################################################
##                                                                                ##
## Initial Preprocessing of Data Including  Incidence Spike Smoothing             ##
##                                                                                ##
##    Ensuring data variables are in correct format, and calculating              ##
##    the incidence of cases outside UK over time.                                ##
##    Also contains code to correct for incidence spike that occurred on          ##
##    May1st. Initial linear interpolation of case numbers between                ##
##    April 30th and 1st May, followed by distribution of cases                   ##
##    equally around 3 days either side.                                          ##
##                                                                                ##
####################################################################################

# Loading in Raw Case Onset Incidence Data
cases_by_onset <- Epicentro_report_onset %>%
  mutate(Date = as.Date(x = Date, format = "%d/%m/%Y")) %>%
  select(Date, Cases_Inside_NorthWest, Cases_UK) %>%
  mutate(Cases_Outside_NorthWest = Cases_UK - Cases_Inside_NorthWest)

# Smoothing of Incidence Spike Around 28th jun
## Creating vector of dates to smooth over
total_days_smoothing <- 2 * days_to_smooth_either_side + 1
dates_vector <- c()
start_date <- smoothing_centre - days_to_smooth_either_side
for (i in 1:total_days_smoothing) {
  dates_vector[i] <- as.character(start_date)
  start_date <- start_date + 1
}

## Linear interpolation of cases on 4th jun from cases on 27th jun & 29th jun
cases_27th_jun <- cases_by_onset[cases_by_onset$Date == "2021-06-27", 2:4]
cases_28th_jun <- cases_by_onset[cases_by_onset$Date == "2021-06-28", 2:4]
cases_29th_jun <- cases_by_onset[cases_by_onset$Date == "2021-06-29", 2:4]
adj_cases_28th_jun <- (cases_27th_jun + cases_29th_jun)/2 # linear intepolation

## Redstributing extra cases to surrounding 3 days either side & 1st May equally 
cases_by_onset[cases_by_onset$Date == "2020-06-28", 2:4] <- adj_cases_28th_jun
cases_to_distribute <- (cases_28th_jun - adj_cases_28th_jun)/total_days_smoothing
rows_for_adjustment <- cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4]
for (i in 1:nrow(rows_for_adjustment)) {
  rows_for_adjustment[i, ] <- rows_for_adjustment[i, ] + cases_to_distribute
}
cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4] <- rows_for_adjustment

####################################################################################
## Calculating Age Distribution of Cases Inside and Outside NorthWest              ##
##                                                                                ##
####################################################################################

# Age Distribution of UK's Population
UK_all_age <- population$proportion[population$location == "Outside"]
UK_all_age <- UK_all_age/100

# Age Distributions of Cases for NorthWest and Nationally 
NorthWest_case_age<- c(8.596, 23.608, 25.277, 17.404, 11.779, 7.838, 3.265, 1.472, 0.761) 
NorthWest_case_age<- NorthWest_case_age/100 # convert percentage to proportion
UK_case_age <- c(6.933, 23.050, 28.661, 17.139, 11.388, 7.700, 3.128, 1.362, 0.644)  
UK_case_age <- UK_case_age/100 # convert percentage to proportion

# Calculating the Age Distribution of Cases Outside NorthWest
NorthWest_Cases_Onset <- cases_by_onset$Cases_Inside_NorthWest # NorthWest only
UK_Cases_Onset <- cases_by_onset$Cases_UK # All of UK, including NorthWest
Outside_NorthWest_Cases <- sum(UK_Cases_Onset) - sum(NorthWest_Cases_Onset)
NorthWest_cases_by_age <- sum(NorthWest_Cases_Onset) * NorthWest_case_age# cases by age for NorthWest
Total_cases_by_age <- sum(UK_Cases_Onset) * UK_case_age # cases by age for UK (inc. NorthWest)
Outside_NorthWest_Cases_by_age <- Total_cases_by_age - NorthWest_cases_by_age # cases by age outside NorthWest
Outside_NorthWest_Case_Age_Dist <- Outside_NorthWest_Cases_by_age/Outside_NorthWest_Cases # proportion cases by age outside NorthWest


####################################################################################
## Calculating Age-Disaggregated Onset Incidence Over Time, for Inside &          ##
##  Outside NorthWest								                                              ##
##                                                                                ##
##    Using the daily onset information for NorthWest and outside NorthWest in    ##
##    conjunction with the age distribution of cases for each of these settings,  ##
##    we are able to calculate onset incidence over time disaggregated by         ##
##    age-group, for each location.                                               ##
##                                                                                ## 
####################################################################################

# Calculating the Incidence of Cases In Different Age Groups Over Time for Inside NorthWest
##  Create dataframe of NorthWest case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_NorthWest <- data.frame(age_groups, NorthWest_case_age)
NorthWest_Age_Cases_Time <- expand.grid(NorthWest_case_age, cases_by_onset$Cases_Inside_NorthWest)
colnames(NorthWest_Age_Cases_Time) <- c("Age_Prop", "Cases")
NorthWest_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(NorthWest_cases_by_age))
NorthWest_Age_Cases_Time <- NorthWest_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_NorthWest, by = c("Age_Prop" = "NorthWest_case_age")) %>%
  mutate(Location = "NorthWest")

# Calculating the Incidence of Cases In Different Age Groups Over Time for Outside NorthWest
##  Create dataframe of outside NorthWest case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_Outside <- data.frame(age_groups, Outside_NorthWest_Case_Age_Dist)
Outside_Age_Cases_Time <- expand.grid(Outside_NorthWest_Case_Age_Dist, cases_by_onset$Cases_Outside_NorthWest)
colnames(Outside_Age_Cases_Time) <- c("Age_Prop", "Cases")
Outside_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(Outside_NorthWest_Cases_by_age))
Outside_Age_Cases_Time <- Outside_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_Outside, by = c("Age_Prop" = "Outside_NorthWest_Case_Age_Dist")) %>%
  mutate(Location = "Outside")

# Combining the Dataframes Together
age_disaggregated_counts_df <- rbind(NorthWest_Age_Cases_Time, Outside_Age_Cases_Time)
age_disaggregated_counts_df <- age_disaggregated_counts_df[, c("Date", "Location", "age_groups", "Cases")]
colnames(age_disaggregated_counts_df) <- c("date", "location", "age_groups", "cases")

# Integrating Age-Disaggregated Cases With Population to Calculate Adjustment Factor
#   Used in downstream analyses
raw_adjustment_factor_df <- age_disaggregated_counts_df %>%
  group_by(age_groups, location) %>%
  summarise(cases = sum(cases)) %>%
  left_join(population, by = c("age_groups" = "age_groups", "location" = "location")) %>%
  mutate(nici = population/cases) %>%
  select(age_groups, location, nici)

# Joining the Case Incidence and Adjustment Factor Datasets Together
age_disaggregated_counts_df <- age_disaggregated_counts_df %>%
  left_join(raw_adjustment_factor_df, by = c("age_groups" = "age_groups", "location" = "location"))

# Saving Created Dataset
saveRDS(age_disaggregated_counts_df, file = "~/desktop/age_disaggregated_onset_incidence_data.rds")

