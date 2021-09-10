# The log likelihood given our parameters

lL_plane <- function(params, x){
  
# Our paramters
  
# mean onset to death
  m_od <- params[1]
# standard deviation of the onet to death
  s_od <- params[2]

  alpha <- 1 / (s_od * s_od)
  beta <- 1 / (m_od * s_od * s_od)
  
# Maximum day of data
  maxday <- params[3]

# Proportion individuals missed in the NorthWest
  z <- params[4]
# Growth rate 
  r <- params[5]
# Detection window
  D <- params[6]
# CFR for each age band
  cfr <- params[7:length(params)]
# adjusting the CFR for each age group, by multiplying each by the highest cfr 
  cfr_adjustment <- cfr[1] * c(1, cfr[2:9])
  cfr <- c(cfr_adjustment[2:9], cfr_adjustment[1])
  
  
# Data unpacking

# Output type
  output <- x[1]
# Number of age bands - 9 in total
  age_bands <- x[2]
# Death observation cut off - varies depending on time period
  death_obv_cut <- x[3]
# "Gold standard" ascertainment - highest level of ascertainment of an age group
  bestz <- x[4]
# Total observed deaths NorthWest
  deaths_x <- x[5]
  deaths_n <- x[6]
# Gov.UK  UK Observed Cases and deaths both total and by age 
  Gov_observed_deaths <- x[7]
  Gov_observed_cases <- x[8]
# Observed deaths by age
  observed_deaths <- x[9:(9 + age_bands - 1)]
  dm <- matrix(x[(9 + age_bands):length(x)], ncol = 5)
# Cases
  cases_subset <- dm[,1]
# Age band
  ageband <- dm[,2]
# Time
  t <- dm[,3]
# Indicator if NorthWest or not
  NorthWest <- dm[,4]
# NiCi
  NiCi <- dm[,5]

  
  
# Model prediction 
# Adjust the number of cases for ascertainment
# This adjustment = 1 for the "gold standard" age group outside of NorthWest
  zNiCi <- ifelse(NorthWest, z, 1 / bestz) * NiCi
  cases_subset_adjusted <- cases_subset * zNiCi
  
# Death gamma density matrix (distributed our death over future days given the O2D)
  dmat <- estimate_deaths_mat(0:71, 72, alpha, beta)
  
# Multiply cases by dmat for (i = NorthWest and outside) and (j = agebands)
  ed <- list()
  k <- 1
  for(i in 0:1){
    for(j in 1:age_bands){
      index <- ageband == j & NorthWest == i
      d2 <- Map("*", cases_subset_adjusted[index] * cfr[j], dmat) 
      ed[[k]] <- Reduce("+", d2)
      k <- k + 1
    }
  }
  pred_deaths <- unlist(ed) # pred deaths is all the deaths that should've occurred in the time window
  pred_deaths[t < death_obv_cut] <- 0 # Adjust for not observing deaths before cutoff
  
# Crude CFR - most basic value
  predicted_deaths_all <- sum(cases_subset_adjusted * cfr[ageband]) # all deaths that should occur based on the cases we've estimated
  observed_cases <- sum(cases_subset)
  crude_cfr <- min(1, predicted_deaths_all/observed_cases)
  
  # Sum cases and deaths by age group
  deaths_by_age <- tapply(pred_deaths, ageband, sum)
  prop_NorthWest <- sum(pred_deaths * NorthWest) / sum(pred_deaths)
  
  ##############################################################################
  
  # Output our predictions if specified
  if(output == 1){
    return(list(
      age_deaths = deaths_by_age,
      prop_NorthWest = prop_NorthWest,
      crude_cfr = crude_cfr
    ))
  }
  
  # Likelihood calculation 
  # Poisson of observed deaths * binomial proportion of deaths in NW * binomial Gov.UK data given the crude CFR
  ll <- sum(dpois(observed_deaths, deaths_by_age, log = TRUE)) +  
    dbinom(deaths_x, deaths_n, prop_NorthWest, log = TRUE) +
    dbinom(Gov_observed_deaths, Gov_observed_cases, crude_cfr, log = TRUE)
  return(ll)
  ##############################################################################
}

# Predict expected distribution of deaths in the future from a given timeperiod and O2D distribution
dp <- function(t, maxday, alpha, beta){
  c(rep(0, t), dgamma(1:(maxday - t), shape = alpha, rate = beta))
}
# Estimate the deaths matrix given the expected distribution of deaths from an O2D distribution 
estimate_deaths_mat <- function(times, maxday, alpha, beta){
  purrr::map(times, dp, maxday = maxday, alpha= alpha, beta = beta)
}


# Prior params from fit to deaths dataset, we spents time on finding the accurate paramter values by 
# using a histogram of the prior data and approximating the curve of best fit over these data points
lP_plane <- function(params){
  m_od <- params[1]
  s_od <- params[2]
  cfr <- params[7]
  
  alpha <- 7
  beta <- 1/0.7
  k <- 10
  shape1 <- 23
  shape2 <- 31
  
  # Prior is gamma, as logical prior of a poisson likelihood 
  if(is.nan(log(beta^alpha/gamma(alpha) * 1/(m_od - k)^(alpha + 1) * exp(-beta/(m_od - k))) )){
    message("m_od", m_od)
    message("s_od", s_od)
  }
  
  # prior is beta * gamma 
  log(beta^alpha/gamma(alpha) * 1/(m_od - k)^(alpha + 1) * exp(-beta/(m_od - k))) +
    dbeta(s_od, shape1, shape2, log = TRUE) +
    8*log(cfr)
  
}

