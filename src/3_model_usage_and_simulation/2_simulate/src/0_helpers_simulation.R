
# 1. Helpers: lag df ------------------------------------------------------

update_covariates_on_day_from_lag_value_df <- function(precip_values_lag, covariates_on_day, day_nr){
  if(!is.null(precip_values_lag$lag3)){
    # 3 days have passed
    # Only inserting if station present and rainfall predicted existed two days ago. Otherwise use existing covariate/lag value from observations.
    if(day_nr == (precip_values_lag$day_nr_lag_1 + 1)){
      covariates_on_day <- covariates_on_day %>% 
        dplyr::left_join(precip_values_lag$lag1, by = c("site" = "site")) %>% 
        dplyr::mutate(lag1_precip = if_else(!is.na(rain_values), rain_values, lag1_precip)) %>% 
        mutate(
          lag1_precip_occurrence = if_else(lag1_precip > 0, 1, 0),
          lag1_precip_shifted = lag1_precip + 0.5,
          log_lag1_precip_shifted = log(lag1_precip_shifted)
        ) %>% dplyr::select(-rain_values)
    }
    if(day_nr == (precip_values_lag$day_nr_lag_2 + 2)){
      covariates_on_day <- covariates_on_day %>% 
        dplyr::left_join(precip_values_lag$lag2, by = c("site" = "site")) %>% 
        dplyr::mutate(lag2_precip = if_else(!is.na(rain_values), rain_values, lag2_precip)) %>% 
        mutate(
          lag2_precip_occurrence = if_else(lag2_precip > 0, 1, 0),
          lag2_precip_shifted = lag2_precip + 0.5,
          log_lag2_precip_shifted = log(lag2_precip_shifted)
        ) %>% dplyr::select(-rain_values)
    }
    if(day_nr == (precip_values_lag$day_nr_lag_3 + 3)){
      covariates_on_day <- covariates_on_day %>% 
        dplyr::left_join(precip_values_lag$lag3, by = c("site" = "site")) %>% 
        dplyr::mutate(lag3_precip = if_else(!is.na(rain_values), rain_values, lag3_precip)) %>% 
        mutate(
          lag3_precip_occurrence = if_else(lag3_precip > 0, 1, 0),
          lag3_precip_shifted = lag3_precip + 0.5,
          log_lag3_precip_shifted = log(lag3_precip_shifted)
        ) %>% dplyr::select(-rain_values)
    }
  }
  
  return(covariates_on_day)
}


update_lag_value_df <- function(precip_values_lag, simulated_rain_values_on_day, covariates_on_day, day_nr){
  
  # Set lag values
  precip_values_lag$lag3 <- precip_values_lag$lag2
  precip_values_lag$day_nr_lag_3 <- precip_values_lag$day_nr_lag_2
  precip_values_lag$lag2 <- precip_values_lag$lag1
  precip_values_lag$day_nr_lag_2 <- precip_values_lag$day_nr_lag_1
  precip_values_lag$lag1 <- data.frame(rain_values = simulated_rain_values_on_day, site = covariates_on_day$site)
  precip_values_lag$day_nr_lag_1 <- day_nr

  return(precip_values_lag)
}


# 2. Mathematical helpers -------------------------------------------------

logit <- function(x){
  return(
    log(x/(1-x))
  )
}

inverse_logit <- function(x){
  return(
    exp(x)/(1+exp(x))
  )
}


# 3. Simulate from gamma distribution -------------------------------------

get_amounts_on_rainy_sites_from_uniform_residuals_GA_gamlss <- function(mu, sigma, uniform_residuals){
  return(
    qGA(
      uniform_residuals,
      mu = mu,
      sigma = sigma
    )
  )
}
