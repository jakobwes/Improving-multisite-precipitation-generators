library(tidyverse)

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
rainy_days <- dat %>% filter(precip > 0)

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
precip ~
  resolution + 
  log_lag1_precip_shifted +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 397063.9, B2-WIC: 396904.5

# Interaction between resolution and lag1:
precip ~
  resolution * 
  log_lag1_precip_shifted +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 397072.5, B2-WIC: 396904.6

# Nonlinearity for lag 1: p spline and sin cos seasonality
precip ~
  resolution + 
  s(log_lag1_precip_shifted) +
  sin_time_of_year +
  cos_time_of_year +
  lat_legendre_polynomial_1 + lat_legendre_polynomial_2 +
  lon_legendre_polynomial_1 + lon_legendre_polynomial_2
# B1_WIC: 397056.4, B2-WIC: 396869.4
# Provisionally taken

# Spatial effect as spline
precip ~
  resolution + 
  s(log_lag1_precip_shifted) +
  sin_time_of_year +
  cos_time_of_year +
  s(lat, lon)
# B1_WIC: 396735.9, B2-WIC: 396571.8
# Taken, improvement in both WIC-definitions

# Seasonality as spline
precip ~
  resolution + 
  s(log_lag1_precip_shifted) +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon)
# B1_WIC: 396651.9, B2-WIC: 396520.4
# Taken, improvement in both WIC-definitions


# Step 1.2 ----------------------------------------------------------------

# Interaction of lag 1 with seasonality:

# Parametric autocorrelation lag 1 -- interaction with parametric seasonality
precip ~
  resolution + 
  log_lag1_precip_shifted*sin_time_of_year + log_lag1_precip_shifted*cos_time_of_year  +
  s(lat, lon)
# B1_WIC: 396584.5, B2-WIC: 396432.3

# Nonparametric autocorrelation lag 1 -- interaction with nonparametric seasonality
precip ~
  resolution + 
  s(log_lag1_precip_shifted, time_of_year_in_days_for_sin_cos)  +
  s(lat, lon)
# B1_WIC: 396424.5, B2-WIC: 396297.6

# Parametric autocorrelation -- interaction with spline seasonality 
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon)
# B1_WIC: 396395.9, B2-WIC: 396305.7
# Taken, lowest WIC in B2 and model much simpler than bivariate spline.

# Lag 2:

# Lag 2 and interaction with seasonality
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon)
# B1_WIC: 395867.6, B2-WIC: 395793.6
# Taken, both WIC lower

# Lag 2 as spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(log_lag2_precip_shifted) +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon)  
# B1_WIC: 396206.8, B2-WIC: 396100.5

# Lag 3 and interaction with seasonality
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon)
# B1_WIC: 395747.2, B2-WIC: 395689.2
# Taken, both WIC lower


# Step 1.3 ----------------------------------------------------------------

precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(lat, lon, time_of_year_in_days_for_sin_cos)
# B1_WIC: 395871.3, B2-WIC: 395590.5
# Not taken, WIC worse


# Step 1.4 ----------------------------------------------------------------

# Altitude: not taken, as it interferes with the estimation of the s(lat, lon) spline and most of the effect is already captured by the lat, lon.

# Alt as linear covariate
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) +
  alt
# B1_WIC: 395755.8, B2-WIC: 395689.5

# Alt in spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) +
  s(alt)
# B1_WIC: 395736.1, B2-WIC: 395682.2
# Taken, WIC improvement


# Step 2 ------------------------------------------------------------------

#   2. Add large scale atmospheric covariates and their internal interactions, as well as their interactions with seasonal effects. Order there
#     2.1. Temperature and interactions with seasonality.
#     2.2. Dewpoint temperature and sea level pressure, as well as their interactions.
#     2.3. Wind speed and wind direction (calculated from the U and V 10m wind component).

# Final model step 1
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) +
  s(alt)
# B1_WIC: 395736.1, B2-WIC: 395682.2


# Step 2.1 ----------------------------------------------------------------

# Temperature as linear factor
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  era5_2m_temperature
# B1_WIC: 395739, B2-WIC: 395676

# Temperature as spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature)
# B1_WIC: 395387.7, B2-WIC: 395310.6
# Taken, lowest B1 and B2 WIC. Interactions with seasonality are probably already captured by the spline

# Temperature as linear factor in interaction with seasonality
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  era5_2m_temperature * sin_time_of_year + era5_2m_temperature * cos_time_of_year
# B1_WIC: 395539.1, B2-WIC: 395450.5

# Temperature as varying coefficient spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(time_of_year_in_days_for_sin_cos, by = era5_2m_temperature, bs = "cc")
# B1_WIC: 743704.6, B2-WIC: 395677


# Step 2.2 ----------------------------------------------------------------

# Mean sea level pressure scaled as linear covariate (using the learnings from parametric model selection that mslp_scaled works better)
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted) +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted) +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted) +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  era5_mean_sea_level_pressure_scaled
# B1_WIC: 387763.3, B2-WIC: 387617.6

# Mean sea level pressure scaled as spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(era5_mean_sea_level_pressure_scaled)
# B1_WIC: 387596.7, B2-WIC: 387537.6

# Mean sea level pressure in interaction with seasonality
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc")
# B1_WIC: 387076.5, B2-WIC: 387024.3
# Chosen, lowest B1 and B2 WIC

# Dewpoint temp difference as linear covariate
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  era5_dewpoint_temp_difference
# B1_WIC: 386507.3, B2-WIC: 386453.3

# Dewpoint temp difference as spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(era5_dewpoint_temp_difference)
# B1_WIC: 386301.1, B2-WIC: 386225.6

# Dewpoint temp difference in interaction with seasonality
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc")
# B1_WIC: 386141.9, B2-WIC: 385994.8
# Chosen, lowest B1 and B2 WIC


# Step 2.3 ----------------------------------------------------------------

# Wind speed as linear effect
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") + 
  era5_wind_speed
# B1_WIC: 386045.3, B2-WIC: 385882.7

# Wind speed as spline
precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
  s(era5_wind_speed)
# B1_WIC: 385816.6, B2-WIC: 385627.5
# Taken, lowest B1 and B2 WIC.


# Final GALM --------------------------------------------------------------
final_galm <- precip ~
  resolution + 
  s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, bs = "cc") +
  s(lat, lon) + s(alt) +
  s(era5_2m_temperature) +
  s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
  s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
  s(era5_wind_speed)

saveRDS(final_galm, "data/4_data_generated/final_models/2_nonparametric/formula_amounts_final_galm.rds")

library("bamlss")
library("gamlss")

fit <- bamlss(
  formula = final_galm,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log"),
  sampler = FALSE
)
saveRDS(fit, "data/4_data_generated/final_models/2_nonparametric/amounts_final_galm.rds")


# Model selection GAMLSS --------------------------------------------------

# Seasonality in sigma
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ cos_time_of_year + sin_time_of_year
)
# B1_WIC: 385427.1, B2-WIC: 385222.3

list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc")
)
# B1_WIC: 385350.9, B2-WIC: 385181.5
# Taken, lower B1 and B2 WIC

# Autocorrelation in sigma with resolution
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + resolution + log_lag1_precip_shifted
)
# B1_WIC: 385137.2, B2-WIC: 384974.7

# Autocorrelation in sigma with resolution and in interaction with seasonality
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc")
)
# B1_WIC: 385092.4, B2-WIC: 384942.3
# Taken, best W1 and W2 BIC

# Autocorrelation in sigma with resolution and in interaction with seasonality 2
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ sin_time_of_year * log_lag1_precip_shifted + cos_time_of_year * log_lag1_precip_shifted + log_lag1_precip_shifted + resolution
)
# B1_WIC: 385223.1, B2-WIC: 384995.7

# Systematic spatial variation
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + s(lat, lon)
)
# B1_WIC: 384962.5, B2-WIC: 384934.1
# Taken, WIC lower

# Wind speed as linear covariate
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    era5_wind_speed
)
# B1_WIC: 384951.5, B2-WIC: 384882.9

# Wind speed as spline
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed)
)
# B1_WIC: 384835.3, B2-WIC: 384811.3
# Taken, lower B1, B2 WIC

# Temperature as linear covariate
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    era5_2m_temperature
)
# B1_WIC: 384859.3, B2-WIC: 384792.1

# Temperature as spline
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature)
)
# B1_WIC: 384698.8, B2-WIC: 384612.1
# Taken, both WIC lower

# Temperature in interaction with seasonality
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(time_of_year_in_days_for_sin_cos, by = era5_2m_temperature, bs = "cc") 
)
# B1_WIC: 392977.8, B2-WIC: 392977.8

# Dewpoint temp difference as linear
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    era5_dewpoint_temp_difference
)
# B1_WIC: 384671.8, B2-WIC: 384585.5

# Dewpoint temp difference as spline
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(era5_dewpoint_temp_difference)
)
# B1_WIC: 384623.9, B2-WIC: 384535.9

# Dewpoint temp difference with seasonality
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc")
)
# B1_WIC: 384584.9, B2-WIC: 384496.7
# Taken lowest WICs

# Mean sea level pressure scaled as linear covariate
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    era5_mean_sea_level_pressure_scaled
)
# B1_WIC: 383982.1, B2-WIC: 383857

# Mean sea level pressure scaled as spline
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_mean_sea_level_pressure_scaled)
)
# B1_WIC: 383978, B2-WIC: 383862
# Not taken because diagnostic plots and simulation indicated substantial model problems

# Mean sea level pressure scaled as linear covariate with seasonality
list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc")
)
# B1_WIC: 383749.8, B2-WIC: 383761.8
# Not taken because diagnostic plots and simulation indicated substantial model problems

# Final BAMLSS ------------------------------------------------------------

final_bamlss <- list(
  precip ~
    resolution + 
    s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag2_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = log_lag3_precip_shifted, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, bs = "cc") +
    s(lat, lon) + s(alt) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_mean_sea_level_pressure_scaled, bs = "cc") +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc") +
    s(era5_wind_speed),
  sigma ~ s(time_of_year_in_days_for_sin_cos, bs = "cc") + 
    resolution + s(time_of_year_in_days_for_sin_cos, by = log_lag1_precip_shifted, bs = "cc") + 
    s(lat, lon) +
    s(era5_wind_speed) +
    s(era5_2m_temperature) +
    s(time_of_year_in_days_for_sin_cos, by = era5_dewpoint_temp_difference, bs = "cc")
)
saveRDS(final_bamlss, "data/4_data_generated/final_models/2_nonparametric/formula_amounts_final_bamlss.rds")

fit <- bamlss(
  formula = final_bamlss,
  data = rainy_days,
  family = GA(mu.link = "log", sigma.link = "log"),
  sampler = FALSE
)
saveRDS(fit, "data/4_data_generated/final_models/2_nonparametric/amounts_final_bamlss.rds")
