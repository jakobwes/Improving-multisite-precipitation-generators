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

effect_types <- "1_parametric"
model_type <- "bamlss"

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

cor_mat <- matrix(1, nr_of_sites, nr_of_sites)
ll_vals <- matrix(1, nr_of_sites, nr_of_sites)
n_datapoints_mat <- matrix(1, nr_of_sites, nr_of_sites)

colnames(cor_mat) <- sites
rownames(cor_mat) <- sites
colnames(ll_vals) <- sites
rownames(ll_vals) <- sites
colnames(n_datapoints_mat) <- sites
rownames(n_datapoints_mat) <- sites

for (i in 1:nr_of_sites) {
  for (j in 1:i) {
    if (i == j) {
      print(paste0("Calculating station ", i))
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
    n_datapoints_mat[i, j] <- n_datapoints
    if (n_datapoints == 0) {
      cor_mat[i, j] <- NA
      ll_vals[i, j] <- NA
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
    
    cor_mat[i, j] <- optim_results$maximum
    ll_vals[i, j] <- optim_results$objective
  }
}

cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
ll_vals[upper.tri(ll_vals)] <- t(ll_vals)[upper.tri(ll_vals)]
n_datapoints_mat[upper.tri(n_datapoints_mat)] <- t(n_datapoints_mat)[upper.tri(n_datapoints_mat)]

saveRDS(cor_mat, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_cor_mat_ml_estimation_latent_field.rds"))
saveRDS(ll_vals, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_ll_vals_ml_estimation_latent_field.rds"))
saveRDS(n_datapoints_mat, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_n_datapoints_mat_ml_estimation_latent_field.rds"))

# 2. Matern infilling --------------------------------------------------------

cor_mat <- readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_cor_mat_ml_estimation_latent_field.rds"))
n_datapoints_mat <- readRDS(paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_n_datapoints_mat_ml_estimation_latent_field.rds"))

sites <- unique(dat$site)
nr_of_sites <- length(sites)
dist_mat <- get_dist_mat(sites, dat)

res_opt <- optim(
  c(0, 0.1, 0),
  cost_function_matern_weighted_all_params,
  cov_mat = cor_mat,
  dist_mat = dist_mat,
  n_datapoints_mat = n_datapoints_mat,
  method = 'BFGS'
)

opt_rho <- exp(res_opt$par[1])
opt_sigma <- res_opt$par[2]
opt_nu <- exp(res_opt$par[3])

cor_mat_matern <- matrix(1, nr_of_sites, nr_of_sites)
colnames(cor_mat_matern) <- sites
rownames(cor_mat_matern) <- sites

for (i in 1:nr_of_sites) {
  for (j in 1:i) {
    if (i == j) {
      cor_mat_matern[i,j] <- 1.0
      next
    }
    
    cov_matern <- get_matern(opt_rho, opt_sigma, dist_mat[i,j], opt_nu)
    
    cor_mat_matern[i,j] <- cov_matern
    cor_mat_matern[j,i] <- cov_matern
  }
}

plot(dist_mat[!diag(nr_of_sites)], cor_mat[!diag(nr_of_sites)], xlab = 'Distance', ylab = 'Correlation', col = rgb(0, 0, 0, alpha = scales::rescale(n_datapoints_mat, to = c(0., 1)))); points(dist_mat[!diag(nr_of_sites)], cor_mat_matern[!diag(nr_of_sites)], col = 'red')

saveRDS(cor_mat_matern, paste0("data/4_data_generated/spatial_dependence/", effect_types, "/", model_type, "_cor_mat_ml_estimation_latent_field_matern.rds"))


# One testing iteration ---------------------------------------------------

site1 <- "S003"
site2 <- "S004"

temp <- get_pivoted_day_classes_for_sites(dat, site1, site2)
days_both_wet <- temp$days_both_wet
days_site1_wet_site2_dry <- temp$days_site1_wet_site2_dry
days_site2_wet_site1_dry <- temp$days_site2_wet_site1_dry
days_both_dry <- temp$days_both_dry
rm(temp)

optimize(
  f = get_log_likelihood,
  interval = c(0.5, 0.99),
  maximum = TRUE,
  days_both_wet = days_both_wet,
  days_site1_wet_site2_dry = days_site1_wet_site2_dry,
  days_site2_wet_site1_dry = days_site2_wet_site1_dry,
  days_both_dry = days_both_dry
)

grid <- seq(0, 0.95, by = 0.025)
values <- sapply(
  grid,
  get_log_likelihood,
  days_both_wet = days_both_wet,
  days_site1_wet_site2_dry = days_site1_wet_site2_dry,
  days_site2_wet_site1_dry = days_site2_wet_site1_dry,
  days_both_dry = days_both_dry
)

max(values)
grid[values == max(values)]

plot(
  grid,
  values
)
