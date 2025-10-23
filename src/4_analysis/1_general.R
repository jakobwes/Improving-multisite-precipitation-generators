rm(list = ls())
library("tidyverse")
library("magrittr")
library("lubridate")
library("reshape2")
library("extRemes")

# SETUP -------------------------------------------------------------------

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
dat %<>% mutate(
  day_nr = as.integer(date - min(date)),
  year = year(date),
  month = month(date)
) %>% arrange(day_nr, site)
rainy_days <- dat %>% filter(precip > 0)


# Map of locations: not the one in the paper due to copyright concerns --------------------------------------------------------

library("ggmap")

locations <- dat %>% distinct(lat, lon, site)

get_map(maptype = "stamen_terrain", bbox = c(bottom = 51, top = 51.5, left = -1.2, right = -0.3), source = "stadia") %>% 
  ggmap(extent = "device") + 
  geom_point(aes(lon, lat, col = "red"), data = locations) + 
  theme(legend.position="none") +
  xlab("Longitude") + ylab("Latitude")

# https://about.google/brand-resource-center/products-and-services/geo-guidelines/
get_map(maptype = "terrain", location = c(lat = 51.25, lon = -0.75)) %>% 
  ggmap(extent = "device") + 
  geom_point(aes(lon, lat, col = "red"), data = locations) + 
  theme(legend.position="none") +
  xlab("Longitude") + ylab("Latitude")
ggsave("locations.png")
