rm(list = ls())
library(gamlss)
library(gamlss.add)
library(tidyverse)
library(magrittr)
library(lubridate)

source("src/2_model_selection/1_parametric_model_selection/0_helpers.R")

# Setup -------------------------------------------------------------------

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
rainy_days <- dat %>% filter(precip > 0)


# Test composite likelihood AIC -------------------------------------------

x <- 1:1000
y <- rnorm(1000, mean = x^2, sd = 10000)
plot(x, y)
wrong_fit <- gamlss(y ~ x)
print_aic_composite_likelihood(wrong_fit)

fit <- gamlss(y ~ poly(x, 2))
print_aic_composite_likelihood(fit)


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
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution + sin_time_of_year + cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,914.3

# Main resolution effect or only interaction:
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted:resolution + log_lag1_precip_shifted + sin_time_of_year + cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,917.4, not taken because main effect seems more principled and AIC higher.

# Add spatial effect: legendre polynomials
# n = 1
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 + lon_legendre_polynomial_1,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,840.3, chosen

# n = 1 + interaction
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,775, chosen

# n = 2
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,745.1, chosen

# n = 2 + interaction
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,744.6, chosen

# n = 3
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 + lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,739.8, chosen.

# n = 3, interaction
fit <- gamlss(
  formula = precip ~ log_lag1_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,708.2, chosen.


# Step 1.2. ---------------------------------------------------------------
#     1.2. Add further autocorrelation effects and interactions between those and seasonality.

# day = 2
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted*resolution + log_lag2_precip_shifted*resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 793,181.2, chosen

# day = 2, day 1 interaction with seasonality
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted*resolution + log_lag2_precip_shifted*resolution + log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,834.9, chosen

# day = 2, day 1 interaction with seasonality with resolution
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted*resolution + log_lag2_precip_shifted*resolution + log_lag1_precip_shifted:sin_time_of_year:resolution + log_lag1_precip_shifted:cos_time_of_year:resolution +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,969.9, not chosen

# day = 2, interaction, day 1 & 2 interaction with seasonality
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,164.4, chosen

# day = 3, interaction
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,127.7, chosen

# day = 3, interaction, interaction with seasonality
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,044.6, chosen


# Step 1.3 ----------------------------------------------------------------

# Altitude
fit <- gamlss(
  formula = precip ~
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    sin_time_of_year + cos_time_of_year +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 792,044.5, chosen


# Step 2. -----------------------------------------------------------------
#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).

# Startpoint: # 792,044.5, chosen

# precip ~ 
    # log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    # log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    # log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    # log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    # sin_time_of_year + cos_time_of_year +
    # lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt

# Temp
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +

    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 791,999.2, but warning from GAMLSS

# Temp and interaction with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 791,796.2, but warning from GAMLSS, chosen

# Mean sea level pressure
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 775,661.8

# Mean sea level pressure with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure + era5_mean_sea_level_pressure:sin_time_of_year + era5_mean_sea_level_pressure:cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 775,467.3, but warning from GAMLSS

# Mean sea level pressure scaled
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 776,285.1, but warning from GAMLSS

# Mean sea level pressure scaled with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 775,403, but warning from GAMLSS. Chosen.

# Dewpoint temperature difference
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 774,291.1, but warning

# Dewpoint temperature difference with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 773,897.4, but warning from GAMLSS. Chosen.

# Wind speed
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 773,661.2, but warning from GAMLSS. Chosen.


# Final GLM  --------------------------------------------------

fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 773,661.2, chosen


# Model selection GAMLSS --------------------------------------------------

# Seasonality in sigma
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year + cos_time_of_year,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,859.7, but warning from GAMLSS, chosen

# Autocorrelation
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year + cos_time_of_year + log_lag1_precip_shifted + resolution,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,686.6

# Autocorrelation occurrence
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year + cos_time_of_year + lag1_precip_occurrence,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,844.3, but warning.

# Autocorrelation with seasonality interaction and resolution
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + resolution,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,663.2, chosen.

# Autocorrelation with seasonality interaction and resolution interaction
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,622.6, chosen.

# Systematic spatial variation: n = 1 + interaction
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,614.6, chosen.

# Systematic spatial variation: n = 2 + n=1 interaction
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 + lon_legendre_polynomial_2,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,611, chosen

# Systematic spatial variation: n = 2 + interactions
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,610.1, chosen

# Wind speed
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~ 
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + 
    era5_wind_speed,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,556, chosen

# Temp
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + 
    era5_wind_speed +
    era5_2m_temperature,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 773,508.1, not chosen

# Temp with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    era5_wind_speed +
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,147.1, but warning, chosen.

# Dewpoint temp difference
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    era5_wind_speed +
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    era5_dewpoint_temp_difference,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,097.8, but warning, chosen

# Dewpoint temp difference with seasonality 
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    era5_wind_speed +
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 772,069.9, but warning, chosen

# Sea level pressure
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    era5_wind_speed +
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    era5_mean_sea_level_pressure_scaled,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 770,899.9, but warning, chosen

# Sea level pressure with seasonality
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    era5_wind_speed +
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)
print_aic_composite_likelihood(fit)
# 770,833.8, but warning, chosen.

# Final model -------------------------------------------------------------
# GLM
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~1,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)


# GAMLSS: AIC chosen model
fit <- gamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  sigma.formula = ~
    # Seasonality and autocorrelation
    sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
    
    # Spatial dependence
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
    
    # Wind speed 
    era5_wind_speed +
    
    # Temperature in interaction with seasonality
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Dewpoint temperature difference in interaction with seasonality
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Mean sea level pressure scaled
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log")
)


# Final model -- BAMLSS -------------------------------------------------------------
library(bamlss)

# For the simulation the same model is used, but fitted by the bamlss-package. 
# It can be verified that the bamlss-fit leads to the same parameter estimates and predictions as the gamlss one, however extracting the variance-covariance matrix 
# (used for the estimation of the AIC_adj) is much more difficult in bamlss. Therefore model selection was done using the gamlss fit.

# BAMLSS -- GLM
fit_bamlss <- bamlss(
  formula = precip ~
    # Autocorrelation and interaction with season
    log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
    log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
    log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
    log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
    
    # Seasonality
    sin_time_of_year + cos_time_of_year +
    
    # Systematic spatial variation
    lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
    
    # Temp
    era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
    
    # Sea level pressure
    era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
    
    # Dewpoint temperature
    era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
    
    # Wind speed
    era5_wind_speed,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log"),
  sampler = FALSE
)
saveRDS(fit_bamlss, "data/4_data_generated/final_models/1_parametric/amounts_final_glm.rds")


# BAMLSS -- GAMLSS (WIC chosen model)
fit_bamlss <- bamlss(
  formula = list(
    precip ~
      # Autocorrelation and interaction with season
      log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
      log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
      log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
      log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
      
      # Seasonality
      sin_time_of_year + cos_time_of_year +
      
      # Systematic spatial variation
      lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
      
      # Temp
      era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
      
      # Sea level pressure
      era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
      
      # Dewpoint temperature
      era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
      
      # Wind speed
      era5_wind_speed,
    sigma ~
      # Seasonality and autocorrelation
      sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
      
      # Spatial dependence
      lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
      
      # Wind speed 
      era5_wind_speed +
      
      # Temperature in interaction with seasonality
      era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
      
      # Dewpoint temperature difference in interaction with seasonality
      era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
      
      # Mean sea level pressure scaled
      era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year
  ),
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log"),
  sampler = FALSE
)
saveRDS(fit_bamlss, "data/4_data_generated/final_models/1_parametric/amounts_with_mslp_final_bamlss.rds")

# Note: when simulating using the model above containing mslp in sigma quite substantial deviations in the tails can be seen. Removing it makes the model a lot better.
# This might be due to the effect form in which it is captured and we did currently not succeed in fitting a parametric form without this tail deviation.
# Therefore we opt to remove it for the final model.
# BAMLSS -- Final model
fit_bamlss <- bamlss(
  formula = list(
    precip ~
      # Autocorrelation and interaction with season
      log_lag1_precip_shifted * resolution + log_lag2_precip_shifted * resolution + log_lag3_precip_shifted*resolution +
      log_lag1_precip_shifted:sin_time_of_year + log_lag1_precip_shifted:cos_time_of_year +
      log_lag2_precip_shifted:sin_time_of_year + log_lag2_precip_shifted:cos_time_of_year +
      log_lag3_precip_shifted:sin_time_of_year + log_lag3_precip_shifted:cos_time_of_year +
      
      # Seasonality
      sin_time_of_year + cos_time_of_year +
      
      # Systematic spatial variation
      lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 + lat_legendre_polynomial_3 * lon_legendre_polynomial_3 + alt +
      
      # Temp
      era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
      
      # Sea level pressure
      era5_mean_sea_level_pressure_scaled + era5_mean_sea_level_pressure_scaled:sin_time_of_year + era5_mean_sea_level_pressure_scaled:cos_time_of_year +
      
      # Dewpoint temperature
      era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year +
      
      # Wind speed
      era5_wind_speed,
    sigma ~
      # Seasonality and autocorrelation
      sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted * resolution +
      
      # Spatial dependence
      lat_legendre_polynomial_1 * lon_legendre_polynomial_1 + lat_legendre_polynomial_2 * lon_legendre_polynomial_2 +
      
      # Wind speed 
      era5_wind_speed +
      
      # Temperature in interaction with seasonality
      era5_2m_temperature + era5_2m_temperature:sin_time_of_year + era5_2m_temperature:cos_time_of_year +
      
      # Dewpoint temperature difference in interaction with seasonality
      era5_dewpoint_temp_difference + era5_dewpoint_temp_difference:sin_time_of_year + era5_dewpoint_temp_difference:cos_time_of_year
  ),
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log"),
  sampler = FALSE
)
saveRDS(fit_bamlss, "data/4_data_generated/final_models/1_parametric/amounts_final_bamlss.rds")



