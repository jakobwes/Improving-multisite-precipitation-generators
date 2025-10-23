
library(mvtnorm)
get_uniform_residuals_standard_multivariate_gaussian <- function(cov_mat){
  return(
    pnorm(
      rmvnorm(
        n = 1,
        mean = rep(0, nrow(cov_mat)),
        sigma = cov_mat
      )
    )
  )
}

get_rainy_site_idxs_from_joint_latent_multivariate_standard_gaussian <- function(p_1, estim_covariance_structure_occurrence, estim_covariance_structure_amounts, estim_covariance_structure_joint, covariates_on_day){
  
  sites_with_observations <- covariates_on_day %>% pull(site)
  simulated_uniform_residuals <- get_uniform_residuals_standard_multivariate_gaussian(as.matrix(estim_covariance_structure_joint[sites_with_observations, sites_with_observations]))
  
  return(
    list(
      "simulated_rainy_sites_idxs_on_day" = which(simulated_uniform_residuals > (1-p_1)),
      "simulated_uniform_residuals" = simulated_uniform_residuals
    )
  )
}

get_rainy_site_idxs_from_seasonal_joint_latent_multivariate_standard_gaussian <- function(p_1, estim_covariance_structure_occurrence, estim_covariance_structure_amounts, estim_covariance_structure_joint, covariates_on_day){
  
  sites_with_observations <- covariates_on_day %>% pull(site)
  season <- covariates_on_day$meteorological_season[1]
  simulated_uniform_residuals <- get_uniform_residuals_standard_multivariate_gaussian(as.matrix(estim_covariance_structure_joint[[season]][sites_with_observations, sites_with_observations]))
  
  return(
    list(
      "simulated_rainy_sites_idxs_on_day" = which(simulated_uniform_residuals > (1-p_1)),
      "simulated_uniform_residuals" = simulated_uniform_residuals
    )
  )
}

get_rain_amounts_on_rainy_sites_from_joint_latent_multivariate_standard_gaussian <- function(predicted_mu, predicted_sigma, estim_covariance_structure_occurrence, estim_covariance_structure_amounts, estim_covariance_structure_joint, covariates_on_day, results_occurrence, p_1){
  simulated_uniform_residuals <- results_occurrence$simulated_uniform_residuals
  simulated_rainy_sites_idxs_on_day <- results_occurrence$simulated_rainy_sites_idxs_on_day
  
  transformed_uniform_residuals_back_to_0_1_for_amounts <- ((simulated_uniform_residuals-(1-p_1))/p_1)[simulated_rainy_sites_idxs_on_day]
  
  return(
    get_amounts_on_rainy_sites_from_uniform_residuals_GA_gamlss(
      mu = predicted_mu, 
      sigma = predicted_sigma,
      uniform_residuals = transformed_uniform_residuals_back_to_0_1_for_amounts
    )
  )
}



