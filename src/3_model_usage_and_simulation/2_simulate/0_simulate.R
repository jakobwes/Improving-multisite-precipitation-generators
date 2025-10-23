library("tidyverse")
library("magrittr")
library("lubridate")
library("reshape2")
library("gamlss")
library("bamlss")

source("src/3_model_usage_and_simulation/2_simulate/src/0_helpers_simulation.R")
source("src/3_model_usage_and_simulation/2_simulate/src/1_joint_occurrence_amounts_structure.R")


simulate <- function(
    dat, 
    occurence_model, 
    amounts_model, 
    
    estim_covariance_structure_occurrence, 
    estim_covariance_structure_amounts, 
    estim_covariance_structure_joint, 
    
    get_rainy_site_idxs_from_positive_probabilities,
    get_simulated_rain_amounts_on_rainy_sites_from_mu_sigma,
    
    verbose = TRUE
){
  
  # Get vector to hold the simulated rain values
  simulated_rain_values <- rep(NA, dim(dat)[1])

  # Data frame in each item with site and simulated precip value
  precip_values_lag <- list(
    lag1 = NULL, 
    day_nr_lag_1 = NULL, 
    lag2 = NULL, 
    day_nr_lag_2 = NULL,
    lag3 = NULL, 
    day_nr_lag_3 = NULL
  )
  
  number_of_days <- length(unique(dat$day_nr))
  i <- 0
  start <- Sys.time()
  
  for(day_nr in unique(dat$day_nr)){
    
    if(i %% 100 == 0 & verbose) print(paste0("Predicting day ", i, " from ", number_of_days, " days, with day_nr ", day_nr))
    
    observation_idxs_on_day <- which(dat$day_nr == day_nr)
    covariates_on_day <- dat[observation_idxs_on_day,]
    simulated_rain_values_on_day <- rep(-1, length(observation_idxs_on_day))
    
    # Update lag values in covariates_on_day. 
    covariates_on_day <- update_covariates_on_day_from_lag_value_df(precip_values_lag, covariates_on_day, day_nr)
    
    # Occurrence: 
    p_1 <- predict(occurrence_model, newdata = covariates_on_day, type = "parameter")
    
    # Amounts: 
    amounts_predictions_on_day <- predict(amounts_model, newdata = covariates_on_day, type = "parameter")
    amounts_predictions_on_day_mu <- amounts_predictions_on_day$mu
    amounts_predictions_on_day_sigma <- amounts_predictions_on_day$sigma
    
    # Get rainy stations from positive probabilities
    results_get_rainy_site_idxs_from_positive_probabilities <- get_rainy_site_idxs_from_positive_probabilities(
      p_1 = p_1,   
      estim_covariance_structure_occurrence = estim_covariance_structure_occurrence, 
      estim_covariance_structure_amounts = estim_covariance_structure_amounts, 
      estim_covariance_structure_joint = estim_covariance_structure_joint,
      covariates_on_day = covariates_on_day
    )
    simulated_rainy_sites_idxs_on_day <- results_get_rainy_site_idxs_from_positive_probabilities$simulated_rainy_sites_idxs_on_day
    
    # Some rainy sites
    if(length(simulated_rainy_sites_idxs_on_day) > 0){
    
      # Set non-rainy sites
      simulated_rain_values_on_day[-simulated_rainy_sites_idxs_on_day] <- 0
      
      # Simulate rainy sites
      simulated_rain_values_on_day[simulated_rainy_sites_idxs_on_day] <- get_simulated_rain_amounts_on_rainy_sites_from_mu_sigma(
        predicted_mu = amounts_predictions_on_day_mu[simulated_rainy_sites_idxs_on_day], 
        predicted_sigma = amounts_predictions_on_day_sigma[simulated_rainy_sites_idxs_on_day],
        estim_covariance_structure_occurrence = estim_covariance_structure_occurrence, 
        estim_covariance_structure_amounts = estim_covariance_structure_amounts, 
        estim_covariance_structure_joint = estim_covariance_structure_joint,
        covariates_on_day = covariates_on_day,
        results_occurrence = results_get_rainy_site_idxs_from_positive_probabilities,
        p_1 = p_1
      )
    } else {
      simulated_rain_values_on_day <- 0
    }
    
    # Save simulated rain values for the day
    simulated_rain_values[observation_idxs_on_day] <- simulated_rain_values_on_day
    
    # Increase index and update lagged precipitation values in df
    i <- i+1
    precip_values_lag <- update_lag_value_df(precip_values_lag, simulated_rain_values_on_day, covariates_on_day, day_nr)

  }

  end <- Sys.time()
  if(verbose) print(paste0("Time elapsed: ", as.numeric(difftime(end, start, units = "mins")), " minutes"))


  return(simulated_rain_values)


}

