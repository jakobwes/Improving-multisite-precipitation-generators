library("tidyverse")
library("magrittr")
library("lubridate")
library("bamlss")
library("gamlss")
library("parallel")

source("src/2_model_selection/2_nonparametric_model_selection/0_bootstrap.R")

# 1. Specify models and arguments -----------------------------------------

# Models to fit and evaluate
formulas <- list(
  precip ~
    resolution + 
    log_lag1_precip_shifted +
    sin_time_of_year +
    cos_time_of_year +
    lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
    lon_legendre_polynomial_1 + lon_legendre_polynomial_2
)

# File to store results
results_filename <- "src/2_model_selection/2_nonparametric_model_selection/results_amounts_1_1.csv"

# Boostrap indices
bootstrapped_dataset_ids <- as.matrix(read.table("src/2_model_selection/2_nonparametric_model_selection/0_bootstrap/bootstrapped_dataset_ids_in_rainy_days.txt"))
n_bootstraps <- if(ncol(bootstrapped_dataset_ids) > 20) 20 else ncol(bootstrapped_dataset_ids)


# 2. Read in data ---------------------------------------------------------

# Data import
dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
rainy_days <- dat %>% filter(precip > 0)


# 3. Helpers --------------------------------------------------------------

# Calculate B1: parameters in the definition of the WIC
get_B1 <- function(vec_ll_bootstrap_model_real_data, vec_ll_boostrap_model_bootstrap_data) {
  # log of / is -
  return(mean(vec_ll_bootstrap_model_real_data - vec_ll_boostrap_model_bootstrap_data))
}

# Calculate B2: parameters in the definition of the WIC
get_B2 <- function(vec_ll_bootstrap_model_real_data, ll_real_model_real_data) {
  # log of / is -
  return(mean(vec_ll_bootstrap_model_real_data - ll_real_model_real_data))
}

# Calculate B1-WIC: negatively oriented so different from Shibata and similar to Ishiguro et al. 1991 we minimize and not maximize
get_WIC_B1 <- function(ll_bootstrap_model_real_data, ll_bootstrap_model_bootstrap_data, ll_real_model_real_data) {
  B1 <- get_B1(ll_bootstrap_model_real_data, ll_bootstrap_model_bootstrap_data)

  return(
    -ll_real_model_real_data - B1
  )
}

# Calculate B2-WIC: negatively oriented so different from Shibata we minimize and not maximize
get_WIC_B2 <- function(ll_bootstrap_model_real_data, ll_real_model_real_data) {
  B2 <- get_B2(ll_bootstrap_model_real_data, ll_real_model_real_data)

  return(
    -ll_real_model_real_data - B2
  )
}

# Calculate the log-likelihoods
get_gamma_ll <- function(values, mu, sigma) {
  return(
    sum(
      log(
        dGA(values, mu = mu, sigma = sigma)
      )
    )
  )
}

# Evaluate on one bootstrap sample
calc_evaluation_on_one_bootstrap_sample <- function(i, formula, bootstrapped_dataset_ids, rainy_days) {
  bootstrapped_model <- bamlss(
    formula = formula,
    family = GA(mu.link = "log", sigma.link = "log"),
    data = rainy_days[bootstrapped_dataset_ids[, i], ],
    sampler = FALSE
  )

  predictions <- predict(bootstrapped_model, newdata = rainy_days, type = "parameter")

  return(list(
    # Model notes
    converged = bootstrapped_model$model.stats$optimizer$converged,
    # LL Bootstrap model, bootstrap data
    ll_boostrap_model_bootstrap_data = logLik(bootstrapped_model)[[1]],
    # LL Bootstrap model real data
    ll_bootstrap_model_real_data = get_gamma_ll(
      values = rainy_days$precip,
      mu = predictions$mu,
      sigma = predictions$sigma
    )
  ))
}
# 5. Evaluate models ------------------------------------------------------

# Fit different GAMLSS-models for mean and scale
n_models <- length(formulas)

models <- tibble(
  model_nr = 1:n_models,
  formula = formulas,
  ll_real_model_real_data = rep(0., n_models),
  training_time = rep(0., n_models),
  bootstraping_time = rep(0., n_models),
  biasedAIC = rep(0., n_models),
  training_converged = FALSE,
  WIC_B1 = rep(0., n_models),
  WIC_B2 = rep(0., n_models)
)

for (temp_i_model in 1:n_models) {
  print(paste0("Training and evaluating model nr. ", temp_i_model))

  formula <- models$formula[temp_i_model][[1]]
  
  # Train real model
  temp_model <- bamlss(
    formula = formula,
    data = rainy_days, 
    family = GA(mu.link = "log", sigma.link = "log"),
    sampler = FALSE
  ) 

  # Log metadata
  models$training_time[temp_i_model] <- as.numeric(temp_model$model.stats$optimizer$runtime)
  models$training_converged[temp_i_model] <- temp_model$model.stats$optimizer$converged

  if (length(AIC(temp_model)) != 0) models$biasedAIC[temp_i_model] <- AIC(temp_model)

  # LL Real model real data
  ll_real_model_real_data <- temp_model$model.stats$optimizer$logPost #as.integer(logLik(temp_model))
  models$ll_real_model_real_data[temp_i_model] <- ll_real_model_real_data

  # Parallelised bootstrap evaluation
  temp_start_bootstraping_time <- Sys.time()

  results_bootstraps <- mclapply(
    X = 1:n_bootstraps,
    FUN = calc_evaluation_on_one_bootstrap_sample,
    formula = formula,
    bootstrapped_dataset_ids = bootstrapped_dataset_ids,
    rainy_days = rainy_days,
    mc.preschedule = TRUE,
    mc.cores = 5
  )

  # Extract results
  bootstrapped_model_converged <- rep(FALSE, n_bootstraps)
  ll_bootstrap_model_bootstrap_data <- rep(0, n_bootstraps)
  ll_bootstrap_model_real_data <- rep(0, n_bootstraps)

  for (temp_j_bootstraping in 1:n_bootstraps) {
    bootstrapped_model_converged[temp_j_bootstraping] <- results_bootstraps[[temp_j_bootstraping]]$converged
    ll_bootstrap_model_bootstrap_data[temp_j_bootstraping] <- results_bootstraps[[temp_j_bootstraping]]$ll_boostrap_model_bootstrap_data
    ll_bootstrap_model_real_data[temp_j_bootstraping] <- results_bootstraps[[temp_j_bootstraping]]$ll_bootstrap_model_real_data
  }

  temp_end_bootstraping_time <- Sys.time()

  # Log bootstraping time as metadata
  models$bootstraping_time[temp_i_model] <- as.numeric(difftime(temp_end_bootstraping_time, temp_start_bootstraping_time, units = "s"))

  # Calculate and save model WICs
  models$WIC_B1[temp_i_model] <- get_WIC_B1(ll_bootstrap_model_real_data, ll_bootstrap_model_bootstrap_data, ll_real_model_real_data)
  models$WIC_B2[temp_i_model] <- get_WIC_B2(ll_bootstrap_model_real_data, ll_real_model_real_data)

  models$formula[temp_i_model] <- paste(deparse(formula), collapse = "")
  
  # Save intermediate results if iteration fails
  write_csv(models, results_filename)
}

# Print results
print(models)


# Save results
write_csv(models, results_filename)
