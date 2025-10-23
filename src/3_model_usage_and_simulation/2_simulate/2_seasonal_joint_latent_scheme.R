library("tidyverse")
library("magrittr")
library("lubridate")
library("reshape2")
library("gamlss")
library("bamlss")

source("src/3_model_usage_and_simulation/2_simulate/0_simulate.R")


# SETUP -------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

effect_types <- args[1] #"2_nonparametric"
model_type <- args[2] #"galm"

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
dat %<>% mutate(day_nr = as.integer(date - min(date)))
dat <- dat %>% arrange(day_nr, site)
rainy_days <- dat %>% filter(precip > 0)

# Start cluster
library(foreach)
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)


# Simulate with joint latent structure --------------------------------------------

# Read in models
amounts_model <- readRDS(paste0("data/4_data_generated/final_models/", effect_types, "/amounts_final_", model_type, ".rds"))
occurrence_model <- readRDS(paste0("data/4_data_generated/final_models/", effect_types, "/occurrence_final_model.rds"))

# Read in joint latent structure
joint_latent_structure <- list(
  "Winter" = readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/seasonal/Winter_", model_type, "_cor_mat_ml_estimation_latent_field_matern.rds")),
  "Spring" = readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/seasonal/Spring_", model_type, "_cor_mat_ml_estimation_latent_field_matern.rds")),
  "Summer" = readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/seasonal/Summer_", model_type, "_cor_mat_ml_estimation_latent_field_matern.rds")),
  "Fall" = readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/seasonal/Fall_", model_type, "_cor_mat_ml_estimation_latent_field_matern.rds"))
)


# Run simulations
simulation_results_joint_latent_scheme <- foreach(1:20, .combine = "cbind", .packages = c("tidyverse", "magrittr", "lubridate", "reshape2", "bamlss", "gamlss", "mvtnorm")) %dopar%{
  
  simulate(
    dat,
    occurence_model = occurrence_model,
    amounts_model = amounts_model,
    
    estim_covariance_structure_occurrence = NULL,
    estim_covariance_structure_amounts = NULL,
    estim_covariance_structure_joint = joint_latent_structure,
    
    get_rainy_site_idxs_from_positive_probabilities = get_rainy_site_idxs_from_seasonal_joint_latent_multivariate_standard_gaussian,
    get_simulated_rain_amounts_on_rainy_sites_from_mu_sigma = get_rain_amounts_on_rainy_sites_from_joint_latent_multivariate_standard_gaussian,
    
    verbose = TRUE
    
  )
}

saveRDS(simulation_results_joint_latent_scheme, paste0("data/4_data_generated/simulation/", effect_types, "/seasonal/", model_type, "_joint_latent_structure.rds"))

# Terminate cluster
stopCluster(cl)
