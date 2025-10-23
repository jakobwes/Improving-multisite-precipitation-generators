rm(list = ls())
library(tidyverse)
library(magrittr)
library(lubridate)
library(reshape2)
library(bamlss)
library(gamlss)

source("src/3_model_usage_and_simulation/1_estimate_spatial_dependence_structure/0_functions.R")

# Setup -------------------------------------------------------------------

output_folder <- "data/4_data_generated/spatial_dependence/"

args <- commandArgs(trailingOnly = TRUE)

effect_types <- args[1] #"2_nonparametric"
model_type <- args[2] #"bamlss"

model_occurrence <- readRDS(paste0("data/4_data_generated/final_models/", effect_types, "/occurrence_final_model.rds"))
model_amounts <- readRDS(paste0("data/4_data_generated/final_models/", effect_types, "/amounts_final_", model_type, ".rds"))

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
dat %<>% mutate(day_nr = as.numeric(date - min(date)))
rainy_days <- dat %>% filter(precip > 0)


# 1. Fit models -----------------------------------------------------------

# Get "residuals"
predictions_occurrence <- predict(model_occurrence, type = "parameter") # if(class(model_occurrence)[1] == "glm") predict(model_occurrence, type = "response") else predict(model_occurrence, type = "parameter")

predictions_amounts <- predict(model_amounts, type = "parameter")
predictions_mu_amounts <- predictions_amounts$mu
predictions_sigma_amounts <- predictions_amounts$sigma

predictions_occurrence_rainy_days <- predictions_occurrence[dat$precip_occurrence]

# Fix params and merge to df
a <- qnorm(1 - predictions_occurrence, 0, 1)
b <- qnorm(
  predictions_occurrence_rainy_days * pGA(rainy_days$precip, mu = predictions_mu_amounts, sigma = predictions_sigma_amounts) + (1 - predictions_occurrence_rainy_days)
)

dat$a <- a
rainy_days$b <- b

dat %<>% left_join(rainy_days %>% dplyr::select(site, date, b), by = c("site", "date"))

# If encountering infinite residual: remove it
if (sum(is.infinite(dat$b)) != 0) dat <- dat[-which(is.infinite(dat$b)), ]
if (sum(is.infinite(rainy_days$b)) != 0) rainy_days <- rainy_days[-which(is.infinite(rainy_days$b)), ]

# Iterate through sites
sites <- unique(dat$site)
nr_of_sites <- length(sites)

cor_mat <- matrix(NA, nr_of_sites, nr_of_sites)
ll_vals <- matrix(NA, nr_of_sites, nr_of_sites)
n_datapoints_mat <- matrix(NA, nr_of_sites, nr_of_sites)

colnames(cor_mat) <- sites
rownames(cor_mat) <- sites
colnames(ll_vals) <- sites
rownames(ll_vals) <- sites
colnames(n_datapoints_mat) <- sites
rownames(n_datapoints_mat) <- sites


# Start cluster
library(foreach)
cl <- parallel::makeCluster(9)
doParallel::registerDoParallel(cl)


res <- simulation_results_joint_latent_scheme <- foreach(i = 1:nr_of_sites, .packages = c("tidyverse", "magrittr", "lubridate", "reshape2", 'mvtnorm')) %dopar%{
  cor_vals <- matrix(nrow = 1, ncol = i)
  ll_vals <- matrix(nrow = 1, ncol = i)
  n_datapoints_vals <- matrix(nrow = 1, ncol = i)
  
  for (j in 1:i) {
    if (i == j) {
      cor_vals[1, j] <- 1.
      ll_vals[1, j] <- NA
      n_datapoints_vals[1, j] <- 1
      next
    }
    
    site1 <- sites[i]
    site2 <- sites[j]
    
    temp <- get_pivoted_day_classes_for_sites(dat, site1, site2)
    days_both_wet <- temp$days_both_wet
    days_site1_wet_site2_dry <- temp$days_site1_wet_site2_dry
    days_site2_wet_site1_dry <- temp$days_site2_wet_site1_dry
    days_both_dry <- temp$days_both_dry
    rm(temp)
    
    n_datapoints <- nrow(days_both_wet) + nrow(days_site1_wet_site2_dry) + nrow(days_site2_wet_site1_dry) + nrow(days_both_dry)
    n_datapoints_vals[1, j] <- n_datapoints
    if (n_datapoints == 0) {
      cor_vals[1, j] <- NA
      ll_vals[1, j] <- NA
      next
    }
    
    optim_results <- optimize(
      f = get_log_likelihood,
      interval = c(0.5, 0.99),
      maximum = TRUE,
      days_both_wet = days_both_wet,
      days_site1_wet_site2_dry = days_site1_wet_site2_dry,
      days_site2_wet_site1_dry = days_site2_wet_site1_dry,
      days_both_dry = days_both_dry
    )
    
    cor_vals[1, j] <- optim_results$maximum
    ll_vals[1, j] <- optim_results$objective
  }
  
  return(list(cor_mat = cor_vals, ll_vals = ll_vals, n_datapoints_vals = n_datapoints_vals))
}
  
# Terminate cluster
stopCluster(cl)

for(i in 1:nr_of_sites){
    cor_mat[i,1:i] <- res[[i]]$cor_mat
    ll_vals[i,1:i] <- res[[i]]$ll_vals
    n_datapoints_mat[i,1:i] <- res[[i]]$n_datapoints_vals
}

cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
ll_vals[upper.tri(ll_vals)] <- t(ll_vals)[upper.tri(ll_vals)]
n_datapoints_mat[upper.tri(n_datapoints_mat)] <- t(n_datapoints_mat)[upper.tri(n_datapoints_mat)]

saveRDS(cor_mat, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_cor_mat_ml_estimation_latent_field.rds"))
saveRDS(ll_vals, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_ll_vals_ml_estimation_latent_field.rds"))
saveRDS(n_datapoints_mat, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_n_datapoints_mat_ml_estimation_latent_field.rds"))
