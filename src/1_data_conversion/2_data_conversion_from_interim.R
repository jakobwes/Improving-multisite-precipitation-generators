rm(list = ls())


# Dependencies ------------------------------------------------------------

library("tidyverse")
library("magrittr")
library("lubridate")


# Load interim data -------------------------------------------------------

dat <- read_csv("data/2_data_interim/interim_data.csv")


# Add shift of precipitation values ---------------------------------------

# All precipitation values below 0.5mm have been already set to zero in the pre-processing.
# We round precipitation values to the nearest 0.1mm. This accounts for the <15 values whose resolution is measured up to second decimal value.
# We fit models to the quantities max(p - 0.5, 0) to account for differences in the capturing in small values.

dat %<>% mutate(
  precip = round(precip, digits = 1),
  precip = if_else(precip - 0.5 > 0, precip - 0.5, 0),
  precip_occurrence = precip > 0
)


# Add covariates ----------------------------------------------------------

# Change in resolution (finer resolution after 1970)
dat %<>% mutate(resolution = if_else(year(date) > 1970, 1, 0))

# Meteorological season
dat %<>% mutate(meteorological_season = cut(dat$month, breaks = c(0, 2, 5, 8, 11, 12), labels = c("Winter", "Spring", "Summer", "Fall", "Winter"), right = TRUE))

# Time of year
dat %<>% mutate(
  time_of_year_in_days = as.numeric(difftime(make_date(year, month, day), make_date(year, 1, 1), units = "days")),
  time_of_year_in_days_for_sin_cos = (time_of_year_in_days / 365) * 2 * pi,
  sin_time_of_year = sin(time_of_year_in_days_for_sin_cos),
  cos_time_of_year = cos(time_of_year_in_days_for_sin_cos)
)

# Lagged values
dat %<>% mutate(date = make_date(year, month, day))

create_lag <- function(df, n_days) {
  column_name <- paste0("lag", n_days, "_precip")
  return(
    df %>% group_by(site) %>% complete(date = full_seq(date, period = 1)) %>% mutate((!!as.symbol(column_name)) := lag(precip, n_days)) %>% ungroup() # Complete if day before not given
  )
}

for (i in 1:4) {
  dat <- create_lag(dat, i)
}
dat %<>% drop_na()

# Lagged values shifted to take the log & lagged occurrence:
lag_columns <- paste0("lag", 1:4, "_precip")
for (column in lag_columns) {
  dat %<>%
    mutate(!!as.symbol(paste0(column, "_occurrence")) := if_else(!!as.symbol(column) > 0, 1, 0)) %>% # TRUE FALSE
    mutate(!!as.symbol(paste0(column, "_shifted")) := !!as.symbol(column) + 0.5)
}

# Log transformation of shifted, lagged values
dat %<>% mutate(across(ends_with("precip_shifted"), ~ log(.x), .names = "log_{col}"))


# Legendre polynomials of lat and lon Source: http://ggorjan.blogspot.com/2009/02/fitting-legendre-orthogonal-polynomials.html
library(orthopolynom)
Legendre <- function(x, n, normalized = TRUE, intercept = FALSE, rescale = TRUE) {
  ## Create a design matrix for Legendre polynomials
  ## x - numeric
  ## n - see orthopolynom
  ## normalized - logical, see orthopolynom
  ## intercept - logical, add intercept
  ## resscale - logical, scale to -1, 1
  if (rescale) x <- scaleX(x, u = -1, v = 1)
  tmp <- legendre.polynomials(n = n, normalized = normalized)
  if (!intercept) tmp <- tmp[2:length(tmp)]
  polynomial.values(polynomials = tmp, x = x, matrix = TRUE)
}

polynomial.values <- function(polynomials, x, matrix = FALSE) {
  ## Changed copy of polynomial.vales from orthopolynom in order
  ## to add matrix argument
  require(polynom)
  n <- length(polynomials)
  if (!matrix) {
    values <- vector(mode = "list", length = n)
  } else {
    values <- matrix(ncol = n, nrow = length(x))
  }
  j <- 1
  while (j <= n) {
    if (!matrix) {
      values[[j]] <- predict(polynomials[[j]], x)
    } else {
      values[, j] <- predict(polynomials[[j]], x)
    }
    j <- j + 1
  }
  values
}

lat_legendre_polynomial <- Legendre(dat$lat, n = 4)
colnames(lat_legendre_polynomial) <- paste0("lat_legendre_polynomial_", 1:4)
lon_legendre_polynomial <- Legendre(dat$lon, n = 4)
colnames(lon_legendre_polynomial) <- paste0("lon_legendre_polynomial_", 1:4)

dat <- cbind(dat, lat_legendre_polynomial, lon_legendre_polynomial)

# Compute some new statistics based on era5

dat %<>% mutate(
  era5_wind_speed = sqrt(era5_10m_u_component_of_wind^2 + era5_10m_v_component_of_wind^2),
  era5_wind_direction = atan(era5_10m_v_component_of_wind / era5_10m_u_component_of_wind),
  era5_dewpoint_temp_difference = era5_2m_temperature-era5_2m_dewpoint_temperature
) %>%
  group_by(month) %>% 
  mutate(era5_mean_sea_level_pressure_scaled = scale(era5_mean_sea_level_pressure)[,1]) %>% 
  ungroup() 

write_csv(dat, "data/3_data_processed/processed_data.csv")

dat %>%
  dplyr::select(
    precip,
    precip_occurrence,
    site,
    date,
    month,
    year,
    alt,
    resolution,
    era5_2m_dewpoint_temperature:era5_mean_sea_level_pressure,
    era5_wind_speed,
    era5_wind_direction,
    era5_dewpoint_temp_difference,
    era5_mean_sea_level_pressure_scaled,
    meteorological_season,
    time_of_year_in_days_for_sin_cos,
    sin_time_of_year,
    cos_time_of_year,
    lag1_precip:lag4_precip,
    lag1_precip_occurrence:lag4_precip_occurrence,
    lag1_precip_shifted:lag4_precip_shifted,
    log_lag1_precip_shifted:log_lag4_precip_shifted,
    lat, 
    lon,
    lat_legendre_polynomial_1:lat_legendre_polynomial_4,
    lon_legendre_polynomial_1:lon_legendre_polynomial_4
  ) %>%
  write_csv("data/3_data_processed/processed_data_covariates_only.csv")
