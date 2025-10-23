library(tidyverse)

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
rainy_days <- dat %>% filter(precip_occurrence > 0)

# Overview:
# Model selection is done using the B1 and B2-WIC following Shibata. We work negatively oriented, so different than from the reference
# we minimize and not maximize. This is consistent with Ishiguro et al. 1991.

# Model selection GLM ------------------------------------------------------

# Steps:
#   1. Fit a base model with lag1 autocorrelation, seasonality and spatial effect
#     1.1. Some selection on spatial effect and seasonality
#     1.2. Add further autocorrelation effects and interactions between those and seasonality.
#     1.3. Add interaction between seasonal and spatial variation.
#     1.4. Add altitude covariates
#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).



# Step 1. -----------------------------------------------------------------

#   1. Fit a base model with lag1 autocorrelation, seasonality and spatial effect
#     1.1. Some selection on spatial effect and seasonality
#     1.2. Add further autocorrelation effects and interactions between those and seasonality.
#     1.3. Add interaction between seasonal and spatial variation.
#     1.4. Add altitude covariates


# Step 1.1 ------------------------------------------------------------------

# Base model with parametric spatial effect and sin cos seasonality:
precip_occurrence ~
  log_lag1_precip_shifted +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 255583.5, B2-WIC: 255572.7

# Nonlinearity for lag 1: p spline and sin cos seasonality
precip_occurrence ~
  s(log_lag1_precip_shifted) +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 254014.6, B2-WIC: 253907.1

# Lag 1 occurrence and amounts as linear effect separately
precip_occurrence ~
  log_lag1_precip_shifted + lag1_precip_occurrence +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 253727.1, B2-WIC: 253674.1

# Lag 1 amounts linear in interaction with seasonality
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 255364, B2-WIC: 255349.9

# Lag 1 amounts linear in interaction with seasonality and occurrence
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + lag1_precip_occurrence +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 253551.4, B2-WIC: 253493.8

# Lag 1 amounts and occurrence linear in interaction with seasonality
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 253464.4, B2-WIC: 253387.1
# Taken because lowest WIC

# Lag 1 amounts as spline and occurrence linear in interaction with seasonality
precip_occurrence ~
  s(log_lag1_precip_shifted) + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 253604.2, B2-WIC: 253504.1

# Lag 1 occurrence linear in interaction with seasonality
precip_occurrence ~
  s(log_lag1_precip_shifted) + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 255073.8, B2-WIC: 255003.7

# Spatial effect as spline
# Lag 1 occurrence and amounts as linear effect separately
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") +
  sin_time_of_year +
  cos_time_of_year +
  s(lat, lon)
# B1_WIC: 253320.5, B2-WIC: 253283.7
# Taken, improvement in both WIC-definitions

# Seasonality as spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 253247.1, B2-WIC: 253229.2
# Taken, improvement in both WIC-definitions


# Step 1.2 ----------------------------------------------------------------

# Autocorrelation lag 2
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted) + lag1_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted) + lag2_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 251429.8, B2-WIC: 251290.7

# Autocorrelation lag 2 occurrence
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 251356.8, B2-WIC: 251305.5

# Autocorrelation lag 2 occurrence and amounts
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 251307.6, B2-WIC: 251216.7
# Taken, improvement in both WIC-definitions

# Autocorrelation lag 3
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted) + lag1_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted) + lag2_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted) + lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 250604.7, B2-WIC: 250382.4

# Autocorrelation lag 3
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + lag3_precip_occurrence+
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 250464.4, B2-WIC: 250318.9

# Autocorrelation lag 3
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 250450.9, B2-WIC: 250314.5
# Taken, improvement in both WIC-definitions

# Autocorrelation lag 3 and interactions (inspired by parametric model)
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 249446.6, B2-WIC: 249294.2

# Autocorrelation lag 3 and interactions (inspired by parametric model)
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 249389.7, B2-WIC: 249236.7

# Autocorrelation lag 3 and interactions (inspired by parametric model)
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon)
# B1_WIC: 249355.6, B2-WIC: 249199.6
# Taken, improvement in both WIC-definitions

# Step 1.3 ----------------------------------------------------------------

precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(lat, lon, time_of_year_in_days_for_sin_cos)
# B1_WIC: 249706.4, B2-WIC: 249120.3
# Not taken, no improvement in WIC

# Step 1.4 ----------------------------------------------------------------

# Alt as linear covariate
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) +
  alt
# B1_WIC: 249363.9, B2-WIC: 249199.8

# Alt in spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) +
  s(alt)
# B1_WIC: 249325.3, B2-WIC: 249199.6


# Step 2 ------------------------------------------------------------------

#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).

# Final model step 1
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) +
  s(alt)
# B1_WIC: 249325.3, B2-WIC: 249199.6


# Step 2.1 ----------------------------------------------------------------

# Temperature as linear factor
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  era5_2m_temperature
# B1_WIC: 262144.3, B2-WIC: 248149.9

# Temperature as spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) +
  s(era5_2m_temperature)
# B1_WIC: 242387.2, B2-WIC: 242113.4
# Taken, lowest B1 WIC. Interactions with seasonality are probably already captured by the spline

# Temperature as varying coefficient spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(time_of_year_in_days_for_sin_cos, by = era5_2m_temperature, bs = "cc")
# B1_WIC: 252583.5, B2-WIC: 249242.6


# Step 2.2 ----------------------------------------------------------------

# Mean sea level pressure scaled as linear effect
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  era5_mean_sea_level_pressure_scaled
# B1_WIC: 218450.5, B2-WIC: 204021.3

# Mean sea level pressure scaled as spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(era5_mean_sea_level_pressure_scaled)
# B1_WIC: 218367.5, B2-WIC: 203898.7

# Mean sea level pressure scaled as varying coefficient spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc")
# B1_WIC: 218286.7, B2-WIC: 203824.1
# Taken, lowest WICs

# Dewpoint temp difference as linear covariate
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  era5_dewpoint_temp_difference
# B1_WIC: 202709.3, B2-WIC: 202049.4

# Dewpoint temp difference as spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(era5_dewpoint_temp_difference)
# B1_WIC: 202344.9, B2-WIC: 201623

# Dewpoint temp difference in interaction with seasonality
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc")
# B1_WIC: 201972.5, B2-WIC: 201274.1
# Taken, best B1 and B2 WIC


# Step 2.3 ----------------------------------------------------------------

# Wind speed as linear effect
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
  era5_wind_speed
# B1_WIC: 200619, B2-WIC: 199962.3

# Wind speed as spline
precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
  s(era5_wind_speed)
# B1_WIC: 200605.5, B2-WIC: 199885.4
# Taken, lowest B1 and B2 WIC.


# Final GALM --------------------------------------------------------------

final_galm <- precip_occurrence ~
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag1_precip_occurrence, bs = "cc") + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag2_precip_occurrence, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") + s(time_of_year_in_days_for_sin_cos, by = lag3_precip_occurrence, bs = "cc") +
  lag1_precip_occurrence : lag2_precip_occurrence + lag2_precip_occurrence : lag3_precip_occurrence + lag1_precip_occurrence : lag2_precip_occurrence : lag3_precip_occurrence +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
  s(era5_wind_speed)

saveRDS(final_galm, "data/4_data_generated/final_models/2_nonparametric/formula_occurrence_final_model.rds")

library("bamlss")

fit <- bamlss(
  formula = final_galm,
  data = dat,
  family = binomial_bamlss,
  sampler = FALSE
)
saveRDS(fit, "data/4_data_generated/final_models/2_nonparametric/occurrence_final_model.rds")
