########################################################################################
##                   September 13th 2020 - November 23rd 2020                         ##
##  1) Data Collation and Preprocessing - An Overview                                 ##
##                                                                                    ##
##   This section loads in data on the demographic structure of the UK population     ##
##   as well as daily case data for NorthWest of the UK. The number of cases are then ##
##   compared against the number of cases in the rest of the UK. We obtain this       ##
##   information from Gov.UK database. This is the UK's government website. This data ##
##   is integrated with information on the Age distribution of the cases. We use      ##
##   this information in order to calculate the age distribution of cases for         ##
##   NorthWest and the rest of UK                                                     ##
##                                                                                    ##
##   Assuming an invariant age-distribution of cases over time, this is then          ##
##   used to convert the aggregate daily incidence of cases (by symptom onset)        ##
##   into age-specific daily case incidences.                                         ##
##                                                                                    ##
########################################################################################

# Loading packages
library(tidyverse)

# Loading in our Data sets from the Gov.UK information
# Some of this data is created from combining datasets
population <- read.csv("~/desktop/UK_population_demography.csv") 
Gov_UK_onset <- read.csv("~/desktop/Gov_UK_data.csv")

# Setting Static Variables of age groups in our 10 year age bands and then smoothing  
# our data at points where extremes occur, we assume these to be outliers
smoothing_centre <- as.Date("2020-10-15", format = "%Y-%m-%d")
days_to_smooth_either_side <- 3
age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

####################################################################################
##                                                                                ##
## Initial Preprocessing of Data Including  Incidence Spike Smoothing             ##
##                                                                                ##
##    Firstly we need to ensure the data is in the correct format                 ##
##    We calculate the incidence of cases outside UK over time.                   ##
##    Also contains code to correct for incidence spike that occurred in          ##
##    early October. Initial linear interpolation of case numbers between         ##
##    October 3rd and October 5th, followed by distribution of cases              ##
##    equally around 3 days either side.                                          ##
##                                                                                ##
####################################################################################

# Loading in Raw Case Onset Incidence Data
cases_by_onset <- Gov_UK_onset %>%
  mutate(Date = as.Date(x = Date, format = "%d/%m/%Y")) %>%
  select(Date, Cases_Inside_NorthWest, Cases_UK) %>%
  mutate(Cases_Outside_NorthWest = Cases_UK - Cases_Inside_NorthWest)

# Smoothing of Incidence Spike Around 04th oct where an outlier occurs, this is the main outlier in all this 
# 72 day period of data
## Creating vector of dates to smooth over
total_days_smoothing <- 2 * days_to_smooth_either_side + 1
dates_vector <- c()
start_date <- smoothing_centre - days_to_smooth_either_side
for (i in 1:total_days_smoothing) {
  dates_vector[i] <- as.character(start_date)
  start_date <- start_date + 1
}

# Linear interpolation of cases on 4th oct from cases on 03th oct & 05th oct
cases_03th_oct <- cases_by_onset[cases_by_onset$Date == "2020-10-03", 2:4]
cases_04th_oct <- cases_by_onset[cases_by_onset$Date == "2020-10-04", 2:4]
cases_05th_oct <- cases_by_onset[cases_by_onset$Date == "2020-10-05", 2:4]

# adjusted the number of cases using the linear interpolation
adj_cases_04th_oct <- (cases_03th_oct + cases_05th_oct)/2 

## Redstributing extra cases to surrounding 3 days either side, as the cases still occured but needed to be evenly 
# distributed as individuals may wait until certain days to go to hospital or be tested
cases_by_onset[cases_by_onset$Date == "2020-10-04", 2:4] <- adj_cases_04th_oct
cases_to_distribute <- (cases_04th_oct - adj_cases_04th_oct)/total_days_smoothing
rows_for_adjustment <- cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4]
for (i in 1:nrow(rows_for_adjustment)) {
  rows_for_adjustment[i, ] <- rows_for_adjustment[i, ] + cases_to_distribute
}
cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4] <- rows_for_adjustment

####################################################################################
## Calculating Age Distribution of Cases Inside and Outside NorthWest             ##
## This requries the use of the proportion  of the cases for both inside and      ##
## and outside of the North West, this is once more obtained from Gov.UK          ##
## databases                                                                      ##
####################################################################################

# Age Distribution of UK's Population
UK_all_age <- population$proportion[population$location == "Outside"]
UK_all_age <- UK_all_age/100

# Age Distributions of Cases for NorthWest and Nationally 
#proportion of cases distributed amongst the North West
NorthWest_case_age<- c(3.950, 12.784, 18.967, 16.880, 14.835, 15.287, 8.458, 4.722, 4.116) 
NorthWest_case_age<- NorthWest_case_age/100 # convert percentage to proportion
#proportion of cases distributed amongst those in the UK, combined from multiple sources data
UK_case_age <- c(4.002, 13.820, 20.174, 16.432, 14.655, 14.745, 7.849, 4.260, 4.064)  
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
##  Calculating Age-Disaggregated Onset of new cases over time, for Inside &      ##
##  Outside NorthWest								                                              ##
##                                                                                ##
##    Using the daily onset information for NorthWest and outside NorthWest in    ##
##    conjunction with the age distribution of cases for each of these settings,  ##
##    we are able to calculate onset incidence over time disaggregated by         ##
##    age-group, for each location.                                               ##
##                                                                                ## 
####################################################################################

# Calculating the Incidence of Cases In Different Age Groups Over Time for Inside NorthWest
# Create dataframe of North West case distribution by age
# Create grid with all combinations of age group and daily case incidence
# Add date to the dataframe
# Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_NorthWest <- data.frame(age_groups, NorthWest_case_age)
NorthWest_Age_Cases_Time <- expand.grid(NorthWest_case_age, cases_by_onset$Cases_Inside_NorthWest)
colnames(NorthWest_Age_Cases_Time) <- c("Age_Prop", "Cases")
NorthWest_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(NorthWest_cases_by_age))
NorthWest_Age_Cases_Time <- NorthWest_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_NorthWest, by = c("Age_Prop" = "NorthWest_case_age")) %>%
  mutate(Location = "NorthWest")

# Calculating the Incidence of Cases In Different Age Groups Over Time for Outside NorthWest
# Create dataframe of outside NorthWest case distribution by age
# Create grid with all combinations of age group and daily case incidence
# Add date to the dataframe
# Multiply cases on each day for by proportion occurring in each age group
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
# Used in downstream analyses
raw_adjustment_factor_df <- age_disaggregated_counts_df %>%
  group_by(age_groups, location) %>%
  summarise(cases = sum(cases)) %>%
  left_join(population, by = c("age_groups" = "age_groups", "location" = "location")) %>%
  mutate(nici = population/cases) %>%
  select(age_groups, location, nici)
#nici is crucial throughout this project as it gives population / cases, the use of this is described in the method section of our report

# Joining the Case Incidence and Adjustment Factor Datasets Together
age_disaggregated_counts_df <- age_disaggregated_counts_df %>%
  left_join(raw_adjustment_factor_df, by = c("age_groups" = "age_groups", "location" = "location"))

# Saving Created Dataset for later use 
saveRDS(age_disaggregated_counts_df, file = "~/desktop/age_disaggregated_onset_incidence_data.rds")
