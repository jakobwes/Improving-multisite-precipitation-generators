rm(list = ls())
library(gamlss)
library(bamlss)
library(tidyverse)
library(magrittr)
library(lubridate)

source("src/2_model_selection/1_parametric_model_selection/0_helpers.R")

# Setup -------------------------------------------------------------------

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")

# Composite likelihood AIC -------------------------------------------

# We use the same adjusted AIC following Varin et al. 2011 for model selection. 
# However vcov(., robust = TRUE) for binomial models returns an error in GAMLSS. So instead a GLM fit is used and 
# model selection is done using estimates from the sandwich package.

# Model selection GLM -----------------------------------------------------

# Steps:
#   1. Fit a base model with lag1 autocorrelation, seasonality and spatial effect
#     1.1. Some selection on spatial effect and levels of polynomials required
#     1.2. Add further autocorrelation effects and interactions between those and seasonality.
#     1.3. Add interaction between seasonal and spatial variation.
#     1.4. Add altitude covariates
#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).


# Step 1.1 ------------------------------------------------------------------

# Base model without spatial effect:
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,248.6

# Add spatial effect: legendre polynomials
# n = 1
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 + lon_legendre_polynomial_1,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,238.6, chosen

# n = 1 + interaction
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,226.7, chosen

# n = 2
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 +
    lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,225.1, chosen

# n = 2 + interaction
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 +
    lat_legendre_polynomial_2 * lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,225.7, not chosen

# n = 3
fit <- glm(
  formula = precip_occurrence ~     
    lag1_precip_occurrence + 
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 +
    lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    lat_legendre_polynomial_1 + lon_legendre_polynomial_3,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 510,226, not chosen

# Step 1.2. ---------------------------------------------------------------
#     1.2. Add further autocorrelation effects and interactions between those and seasonality.

# day = 1 amounts with resolution as linear effect
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 508,298.3, chosen

# day = 1 amounts interaction with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    lag1_precip_occurrence + 
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 507,926.8, chosen

# day = 1 amounts and occurence interaction with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 507,836.7, chosen

# day = 2
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    lag2_precip_occurrence + lag2_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 503,718.4, chosen

# day = 2, interaction
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 503,604.1, chosen

# day = 2, interaction of occurrences
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 501,453.2, chosen

# day = 3,
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 499,744.3, chosen

# day = 3,
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + 
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 499,632.1, chosen

# day = 3,
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 499,551.3, chosen



# Step 1.3. ---------------------------------------------------------------
#     1.3. Add interaction between seasonal and spatial variation.

# Interaction with degree 1 polynomial
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_1 * sin_time_of_year + lon_legendre_polynomial_1 * sin_time_of_year+ lat_legendre_polynomial_2 * sin_time_of_year + lon_legendre_polynomial_2 * sin_time_of_year +
    lat_legendre_polynomial_1 * cos_time_of_year + lon_legendre_polynomial_1 * cos_time_of_year+ lat_legendre_polynomial_2 * cos_time_of_year + lon_legendre_polynomial_2 * cos_time_of_year,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 499,553.2 not pursued further, because it adds a lot of model complexity without gains


# Step 2. -----------------------------------------------------------------
#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).

# Startpoint: # 499,551.3, chosen

# fit <- glm(
#   formula = precip_occurrence ~ 
      # sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
      # sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
      # sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
      # lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
      # sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
      # sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
      # sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
      # lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
#   data = dat,
#   family = binomial
# )


# Temp
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    era5_2m_temperature,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 497,474, chosen

# Temp and interaction with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 491,412.1, chosen

# Mean sea level pressure 
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 417,275.4

# Mean sea level pressure in interaction with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure*sin_time_of_year + era5_mean_sea_level_pressure*cos_time_of_year,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 413,624.2

# Mean sea level pressure scaled
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 412,638.1

# Mean sea level pressure scaled in interaction with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 412,521.5, chosen

# Dewpoint temperature
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    era5_2m_dewpoint_temperature,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 408,276.1

# Dewpoint temperature difference
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    era5_dewpoint_temp_difference,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 408,276.1

# Dewpoint temperature difference with seasonality
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    era5_dewpoint_temp_difference*sin_time_of_year + era5_dewpoint_temp_difference*cos_time_of_year,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 407,500, chosen

# Wind speed
fit <- glm(
  formula = precip_occurrence ~ 
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    era5_dewpoint_temp_difference*sin_time_of_year + era5_dewpoint_temp_difference*cos_time_of_year +
    era5_wind_speed,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# 404,079.1, chosen


# Final model -------------------------------------------------------------

fit <- glm(
  formula = precip_occurrence ~ 
    # Autocorrelation occurrence
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    
    # Autocorrelation amounts
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    
    # Spatial dependence
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    
    # Temperature
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    
    # Mean sea level pressure scaled
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference*sin_time_of_year + era5_dewpoint_temp_difference*cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  data = dat,
  family = binomial
)
print_sandwich_aic_adj(fit)
# AIC_adj = 404,079.1

# Final model in bamlss (for computational consistency)
fit <- bamlss(
  precip_occurrence ~ 
    # Autocorrelation occurrence
    sin_time_of_year*lag1_precip_occurrence + cos_time_of_year*lag1_precip_occurrence +
    sin_time_of_year*lag2_precip_occurrence + cos_time_of_year*lag2_precip_occurrence +
    sin_time_of_year*lag3_precip_occurrence + cos_time_of_year*lag3_precip_occurrence +
    lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
    
    # Autocorrelation amounts
    sin_time_of_year*lag1_precip_shifted + cos_time_of_year*lag1_precip_shifted + lag1_precip_shifted*resolution +
    sin_time_of_year*lag2_precip_shifted + cos_time_of_year*lag2_precip_shifted + lag2_precip_shifted*resolution +
    sin_time_of_year*lag3_precip_shifted + cos_time_of_year*lag3_precip_shifted + lag3_precip_shifted*resolution +
    
    # Spatial dependence
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2 +
    
    # Temperature
    sin_time_of_year*era5_2m_temperature + 
    cos_time_of_year*era5_2m_temperature +
    
    # Mean sea level pressure scaled
    era5_mean_sea_level_pressure_scaled*sin_time_of_year + era5_mean_sea_level_pressure_scaled*cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference*sin_time_of_year + era5_dewpoint_temp_difference*cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  data = dat,
  family = binomial_bamlss,
  sampler = FALSE
)
saveRDS(fit, "data/4_data_generated/final_models/1_parametric/occurrence_final_model.rds")



