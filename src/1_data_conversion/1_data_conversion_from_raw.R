rm(list = ls())

# Dependencies ------------------------------------------------------------

library("Rglimclim")
library("tidyverse")
library("magrittr")
library("lubridate")
library(ncdf4) # needed to read netcdf data
library(ncdf4.helpers) # additional support functions for netcdf data


# Read in Blackwater dataset ----------------------------------------------

load("data/1_data_raw/BlackwaterData/BlackwaterSiteinfo.rda")
dat <- as_tibble(read.GLCdata("data/1_data_raw/BlackwaterData/BlackwaterData_SelectedV2.dat"))


# Reformat data and merge siteinfo ----------------------------------------

dat %<>% rename(
  year = Year,
  month = Month,
  day = Day,
  site = Site,
  precip = Var1
)

site.frame %<>% rename(
  lat = Latitude,
  lon = Longitude,
  alt = Altitude
)

dat %<>% left_join(site.frame %>% dplyr::select(SCode, lat, lon, alt), by = c("site" = "SCode"))

# Remove combinations without metadata:
dat %<>% drop_na()


# Add ERA5 covariates -----------------------------------------------------
# c("2m_temperature", "2m_dewpoint_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "mean_sea_level_pressure")

folder <- "data/1_data_raw/MeteorologicalVariables/"

# 2m_dewpoint_temperature
name <- "era5_2m_dewpoint_temperature"

nc_object <- nc_open(paste0(folder, name, ".nc"))
variable <- ncvar_get(nc_object, "d2m")
time <- as.character(nc.get.time.series(nc_object))
nc_close(nc_object)
remove(nc_object)

era5_data <- data.frame("variable" = variable, "date" = time) %>% rename(!!as.symbol(name) := variable)


# 2m_temperature
name <- "era5_2m_temperature"

nc_object <- nc_open(paste0(folder, name, ".nc"))
variable <- ncvar_get(nc_object, "t2m")
time <- as.character(nc.get.time.series(nc_object))
nc_close(nc_object)
remove(nc_object)

era5_data %<>% left_join(
  data.frame("variable" = variable, "date" = time),
  by = "date"
) %>% rename(!!as.symbol(name) := variable)


# 10m_u_component_of_wind
name <- "era5_10m_u_component_of_wind"

nc_object <- nc_open(paste0(folder, name, ".nc"))
variable <- ncvar_get(nc_object, "u10")
time <- as.character(nc.get.time.series(nc_object))
nc_close(nc_object)
remove(nc_object)

era5_data %<>% left_join(
  data.frame("variable" = variable, "date" = time),
  by = "date"
) %>% rename(!!as.symbol(name) := variable)


# 10m_v_component_of_wind
name <- "era5_10m_v_component_of_wind"

nc_object <- nc_open(paste0(folder, name, ".nc"))
variable <- ncvar_get(nc_object, "v10")
time <- as.character(nc.get.time.series(nc_object))
nc_close(nc_object)
remove(nc_object)

era5_data %<>% left_join(
  data.frame("variable" = variable, "date" = time),
  by = "date"
) %>% rename(!!as.symbol(name) := variable)


# mean_sea_level_pressure
name <- "era5_mean_sea_level_pressure"

nc_object <- nc_open(paste0(folder, name, ".nc"))
variable <- ncvar_get(nc_object, "msl")
time <- as.character(nc.get.time.series(nc_object))
nc_close(nc_object)
remove(nc_object)

era5_data %<>% left_join(
  data.frame("variable" = variable, "date" = time),
  by = "date"
) %>% rename(!!as.symbol(name) := variable)


# Do some transformations and merge onto df
era5_data %<>% mutate(
  date = as.Date(date)
) %>% mutate(day = day(date), month = month(date), year = year(date))

dat %<>% left_join(era5_data, by = c("day", "month", "year"))

# Remove combinations without era5 data
dat %<>% drop_na()


# Remove S089 as there are data issues there ------------------------------

dat %<>% filter(site != "S089")


# Write to file -----------------------------------------------------------

write_csv(dat, "data/2_data_interim/interim_data.csv")
