rm(list = ls())
library("tidyverse")
library("magrittr")
library("lubridate")
library("reshape2")
library("extRemes")

# SETUP -------------------------------------------------------------------

effect_types <- "1_parametric"
model_type_loc <- if (effect_types == "2_nonparametric") "galm" else "glm"
model_name_loc <- if (effect_types == "2_nonparametric") "GAM" else "GLM"

dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv")
dat %<>% mutate(
  day_nr = as.integer(date - min(date)),
  year = year(date),
  month = month(date)
) %>% arrange(day_nr, site)
rainy_days <- dat %>% filter(precip > 0)

results_galm <- data.frame(readRDS(file = paste0("data/4_data_generated/simulation/", effect_types, "/seasonal/", model_type_loc, "_joint_latent_structure.rds")))
results_bamlss <- data.frame(readRDS(file = paste0("data/4_data_generated/simulation/", effect_types, "/seasonal/bamlss_joint_latent_structure.rds")))

colnames(results_galm) <- paste0("galm_simulation_run_", c(1:ncol(results_galm)))
colnames(results_bamlss) <- paste0("bamlss_simulation_run_", c(1:ncol(results_bamlss)))

# Add enriching information
results_galm <- results_galm %>% mutate(
  date = dat$date,
  site = dat$site,
  meteorological_season = dat$meteorological_season,
  real_precip = dat$precip
)

results_bamlss <- results_bamlss %>% mutate(
  date = dat$date,
  month = dat$month,
  year = dat$year,
  site = dat$site,
  meteorological_season = dat$meteorological_season,
  real_precip = dat$precip
)

# Pivot for calculation
results_galm <- results_galm %>% pivot_longer(cols = starts_with("galm_simulation_run_"), names_to = "simulation_run", values_to = "simulation_values")
results_bamlss <- results_bamlss %>% pivot_longer(cols = starts_with("bamlss_simulation_run_"), names_to = "simulation_run", values_to = "simulation_values")

# Full merge
results <- cbind(results_galm %>% rename(simulation_run_galm = simulation_run, simulation_values_galm = simulation_values), results_bamlss %>% rename(simulation_run_bamlss = simulation_run, simulation_values_bamlss = simulation_values))
results <- results[, unique(colnames(results))]

# Standardize columns
results %<>% mutate(
  simulation_run = substr(simulation_run_galm, 21, nchar(simulation_run_galm))
) %>%
  dplyr::select(-c(simulation_run_galm, simulation_run_bamlss))

# Subset to 19 simulation runs
results %<>% filter(simulation_run != 20)


# Helpers  ----------------------------------------------------------------


make_qq_plot_of_results <- function(results, ymax = 120, main = "", model_type_loc = "GALM"){

  model_type_loc_name <- toupper(model_type_loc)

  results_qq_plot <- results %>% 
    group_by(simulation_run) %>% 
    mutate(
      simulation_values_galm = sort(simulation_values_galm),
      simulation_values_bamlss = sort(simulation_values_bamlss),
      real_precip = sort(real_precip),
      id = 1:length(real_precip)
    ) %>% 
    ungroup() %>% group_by(id) %>% 
    summarise(
      real_precip = first(real_precip),
      galm_med = median(simulation_values_galm), bamlss_med = median(simulation_values_bamlss), 
      galm_min = min(simulation_values_galm), bamlss_min = min(simulation_values_bamlss), 
      galm_max = max(simulation_values_galm), bamlss_max = max(simulation_values_bamlss)
    ) %>% ungroup() 
  
  results_qq_plot %>% ggplot() +
    coord_cartesian(ylim = c(0, ymax), xlim = c(0, ymax)) +
    xlab("Observed Quantiles (in mm)") + ylab("Simulated Quantiles (in mm)") +
    scale_color_manual(
      values = setNames(c("indianred1", "lightblue"), c(model_type_loc_name, "GAMLSS")),
      breaks = c(model_type_loc_name, "GAMLSS")
    ) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120)) +
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120)) +
    ggtitle(main) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "orange", linetype = 2) +
    geom_pointrange(aes(x = real_precip, y = galm_med, ymin = galm_min, ymax = galm_max, colour = model_type_loc_name), alpha = 0.6) + 
    geom_pointrange(aes(x = real_precip, y = bamlss_med, ymin = bamlss_min, ymax = bamlss_max, colour = "GAMLSS"), alpha = 0.6) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.title=element_blank(),
                       axis.title = element_text(size = 16),   # axis labels
                       axis.text = element_text(size = 14),    # tick labels
                       legend.text = element_text(size = 16),  # legend text
    ) +guides(colour = guide_legend(nrow = 1))
}


# QQ-Plots: Values ----------------------------------------------------------------

# Overall
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/"))
make_qq_plot_of_results(results, model_type_loc = model_name_loc)
ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/overall.png"), width = 1800, height = 1400, units = "px")


# Seasonal
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/meteorological_season/"))
for (it_meteorological_season in unique(results$meteorological_season)) {
  
  make_qq_plot_of_results(results %>% filter(meteorological_season == it_meteorological_season), model_type_loc = model_name_loc)
  ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/meteorological_season/", it_meteorological_season, ".png"), width = 1800, height = 1400, units = "px")
  
}


# Monthly
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/months/"))
month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
for (it_month in unique(results$month)){
  
  make_qq_plot_of_results(results %>% filter(month == it_month), model_type_loc = model_name_loc)
  ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/months/", it_month, ".png"), width = 1800, height = 1400, units = "px")
  
}


# Sitewide
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites/"))
for (it_site in unique(results$site)) {
  
  make_qq_plot_of_results(results %>% filter(site == it_site), model_type_loc = model_name_loc)
  ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites/", it_site, ".png"), width = 1800, height = 1400, units = "px")
  
}


# Seasonal, sitewide
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_seasons/"))
for (it_meteorological_season in unique(results$meteorological_season)) {
  dir.create(file.path("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_seasons", it_meteorological_season), showWarnings = FALSE)
  for (it_site in unique(results$site)) {
    
    make_qq_plot_of_results(results %>% filter(site == it_site, meteorological_season == it_meteorological_season), model_type_loc = model_name_loc)
    ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_seasons/", it_meteorological_season, "/", it_site, ".png"), width = 1800, height = 1400, units = "px")
    
  }
}


# Monthly, sitewide
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_months/"))
month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
for (it_month in 1:12) {
  dir.create(file.path("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_months", month_names[it_month]), showWarnings = FALSE)
  for (it_site in unique(results$site)) {
    
    make_qq_plot_of_results(results %>% filter(site == it_site, month == it_month), model_type_loc = model_name_loc)
    ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/sites_months/", month_names[it_month], "/", it_site, ".png"), width = 1800, height = 1400, units = "px")
    
  }
}


# QQ-plots: Mean rainfall ----------------------------------------------

## Get mean rainfall values
mean_precip_df <- results %>% 
  group_by(simulation_run, date, month, meteorological_season) %>% 
  summarise(
    real_precip = mean(real_precip),
    simulation_values_galm = mean(simulation_values_galm),
    simulation_values_bamlss = mean(simulation_values_bamlss),
    .groups = "drop"
  )



# Overall
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/"))
make_qq_plot_of_results(mean_precip_df, model_type_loc = model_name_loc)
ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/overall.png"), width = 1800, height = 1400, units = "px")


# Seasonal
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/meteorological_season/"))
for (it_meteorological_season in unique(mean_precip_df$meteorological_season)) {
  
  make_qq_plot_of_results(mean_precip_df %>% filter(meteorological_season == it_meteorological_season), model_type_loc = model_name_loc)
  ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/meteorological_season/", it_meteorological_season, ".png"), width = 1800, height = 1400, units = "px")
  
}


# Monthly
dir.create(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/months/"))
month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
for (it_month in unique(mean_precip_df$month)){
  
  make_qq_plot_of_results(mean_precip_df %>% filter(month == it_month), model_type_loc = model_name_loc)
  ggsave(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/mean_rainfall/months/", it_month, ".png"), width = 1800, height = 1400, units = "px")
  
}



# Spatial dependence Occurrence ------------------------------------------------------

# Histogram percent rainy stations

# Observations
png(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/hist_occurrence_obs.png"), width = 600, height = 500)
results %>% 
  filter(simulation_run == 1) %>%
  group_by(date) %>%
  summarise(
    quotient = sum((real_precip>0))/dplyr::n(), .groups = "drop"
  ) %>% 
  pull(quotient) %>% 
  hist(freq = FALSE, breaks = 10, main = "Observations")
dev.off()

# GALM
png(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/hist_occurrence_galm.png"), width = 600, height = 500)
results %>% 
  group_by(date, simulation_run) %>%
  summarise(
    quotient = sum((simulation_values_galm>0))/dplyr::n(), .groups = "drop"
  ) %>% 
  pull(quotient) %>% 
  hist(freq = FALSE, breaks = 10, main = toupper(model_type_loc))
dev.off()

# GAMLSS
png(paste0("plots/", effect_types, "/joint_latent_multivariate_gaussian/seasonal/hist_occurrence_bamlss.png"), width = 600, height = 500)
results %>% 
  group_by(date, simulation_run) %>%
  summarise(
    quotient = sum((simulation_values_bamlss>0))/dplyr::n(), .groups = "drop"
  ) %>% 
  pull(quotient) %>% 
  hist(freq = FALSE, breaks = 10, main = "GAMLSS")
dev.off()


# Spatial dependence amounts --------------------------------------------------

# Amounts correlation in real data
precip_matrix <- reshape2::dcast(rainy_days, date ~ site, value.var = "precip")
precip_matrix <- precip_matrix[, - which(colnames(precip_matrix) == "date")]
correlation <- cor(precip_matrix, use = "pairwise.complete.obs")
constant_intersite_correlation_amounts <- mean(correlation, na.rm = TRUE)

# Bamlss-Amounts correlation in one simulation run
precip_matrix_bamlss <- reshape2::dcast(results_bamlss %>% filter(simulation_values > 0, simulation_run == "bamlss_simulation_run_4"), date ~ site, value.var = "simulation_values")
precip_matrix_bamlss <- precip_matrix_bamlss[, - which(colnames(precip_matrix_bamlss) == "date")]
correlation_bamlss <- cor(precip_matrix_bamlss, use = "pairwise.complete.obs")
constant_intersite_correlation_bamlss <- mean(correlation_bamlss, na.rm = TRUE)

# Print results of comparison with one run
print(mean((correlation_bamlss - correlation)^2, na.rm = TRUE))
print(constant_intersite_correlation_bamlss)
print(constant_intersite_correlation_amounts)

# galm-Amonts correlation in one simulation run
precip_matrix_galm <- reshape2::dcast(results_galm %>% filter(simulation_values > 0, simulation_run == "galm_simulation_run_4"), date ~ site, value.var = "simulation_values")
precip_matrix_galm <- precip_matrix_galm[, - which(colnames(precip_matrix_galm) == "date")]
correlation_galm <- cor(precip_matrix_galm, use = "pairwise.complete.obs")
constant_intersite_correlation_galm <- mean(correlation_galm, na.rm = TRUE)

# Print results of comparison with one run
print(mean((correlation_galm - correlation)^2, na.rm = TRUE))
print(constant_intersite_correlation_galm)
print(constant_intersite_correlation_amounts)


# Correlation in bamlss-Amounts as histogram and distribution across runs:
# Real data
precip_matrix <- reshape2::dcast(rainy_days, date ~ site, value.var = "precip")
precip_matrix <- precip_matrix[, - which(colnames(precip_matrix) == "date")]
correlation <- cor(precip_matrix, use = "pairwise.complete.obs")
constant_intersite_correlation_amounts <- mean(correlation, na.rm = TRUE)

# bamlss values
cor_vals <- rep(0, 19)
mean_squared_distance <- rep(0, 19)
for(i in 1:19){
  precip_matrix_bamlss <- reshape2::dcast(results_bamlss %>% filter(simulation_values > 0, simulation_run == paste0("bamlss_simulation_run_",i)), date ~ site, value.var = "simulation_values")
  precip_matrix_bamlss <- precip_matrix_bamlss[, - which(colnames(precip_matrix_bamlss) == "date")]
  correlation_bamlss <- cor(precip_matrix_bamlss, use = "pairwise.complete.obs")
  mean_squared_distance[i] <- mean((correlation_bamlss - correlation)^2, na.rm = T)
  cor_vals[i] <- mean(correlation_bamlss, na.rm = TRUE)
}

print(mean(cor_vals))
print(sd(cor_vals))
print(mean(mean_squared_distance))

# Season-amounts-bamlss correlation across runs
constant_intersite_correlation_amounts <-  list(Spring = 0, Summer = 0, Fall = 0, Winter = 0)

cor_vals <- data.frame(Spring = rep(0, 19), Summer = rep(0, 19), Fall = rep(0, 19), Winter = rep(0, 19))
mean_squared_distance <- data.frame(Spring = rep(0, 19), Summer = rep(0, 19), Fall = rep(0, 19), Winter = rep(0, 19))

for(season in colnames(cor_vals)){
  # Real data seasonal correlations
  precip_matrix <- reshape2::dcast(rainy_days %>% filter(meteorological_season == season), date ~ site, value.var = "precip")
  precip_matrix <- precip_matrix[, - which(colnames(precip_matrix) == "date")]
  correlation <- cor(precip_matrix, use = "pairwise.complete.obs")
  constant_intersite_correlation_amounts[season] <- mean(correlation, na.rm = TRUE)
  
  # Bamlss simulations seasonal correlations
  for(i in 1:19){
    precip_matrix_bamlss <- reshape2::dcast(results_bamlss %>% filter(simulation_values > 0, meteorological_season == season, simulation_run == paste0("bamlss_simulation_run_",i)), date ~ site, value.var = "simulation_values")
    precip_matrix_bamlss <- precip_matrix_bamlss[, - which(colnames(precip_matrix_bamlss) == "date")]
    correlation_bamlss <- cor(precip_matrix_bamlss, use = "pairwise.complete.obs")
    mean_squared_distance[i, season] <- mean((correlation_bamlss - correlation)^2, na.rm = T)
    cor_vals[i, season] <- mean(correlation_bamlss, na.rm = TRUE)
  }
}

apply(cor_vals, 2, mean)
constant_intersite_correlation_amounts
apply(mean_squared_distance, 2, mean)


# Occurrence analysis -----------------------------------------------------

results %>% 
  group_by(simulation_run) %>% 
  summarise(occ_real = mean(real_precip > 0), occ_galm = mean(simulation_values_galm > 0), occ_bamlss = mean(simulation_values_bamlss > 0), .groups = "drop") %>% 
  summarise(occ_real = mean(occ_real), sd_occ_galm = sd(occ_galm), sd_occ_bamlss = sd(occ_bamlss), occ_galm = mean(occ_galm), occ_bamlss = mean(occ_bamlss))

results %>% 
  group_by(simulation_run, site) %>% 
  summarise(occ_real = mean(real_precip > 0), occ_galm = mean(simulation_values_galm > 0), occ_bamlss = mean(simulation_values_bamlss > 0), .groups = "drop") %>% 
  mutate(mae_galm = abs(occ_real - occ_galm), mae_bamlss = abs(occ_real - occ_bamlss)) %>%
  summarise(sd_mae_galm = sd(mae_galm), mae_galm = mean(mae_galm), sd_mae_bamlss = sd(mae_bamlss), mae_bamlss = mean(mae_bamlss))


# GEV Analysis -- Yearly maxima --------------------------------------------------------------

sites_for_gev <- results %>% group_by(site) %>% summarise(n = n_distinct(year), .groups= "drop") %>% filter(n >= 30) %>% pull(site)

yearly_max_by_site <- results %>%
  filter(site %in% sites_for_gev) %>%
  group_by(
    year,
    site,
    simulation_run
  ) %>%
  summarise(
    yearly_max_real = max(real_precip),
    yearly_max_galm = max(simulation_values_galm),
    yearly_max_bamlss = max(simulation_values_bamlss),
    .groups = "drop"
  )

yearly_max_by_site <- rbind(yearly_max_by_site, yearly_max_by_site %>% mutate(site = "pooled"))

library("extRemes")

vec_qevd <- Vectorize(qevd)

fit_fevd_and_get_params <- function(data) {
  fit <- fevd(data, type = "GEV", method = "MLE")
  return(tibble(
    location = distill.fevd(fit)["location"],
    scale = distill.fevd(fit)["scale"],
    shape = distill.fevd(fit)["shape"]
  ))
}

get_quantiles_from_params_of_data_type <- function(.data, type) {
  mutate(
    .data,
    !!as.symbol(paste0("q50_", type)) := vec_qevd(0.5, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV"),
    !!as.symbol(paste0("q90_", type)) := vec_qevd(0.9, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV"),
    !!as.symbol(paste0("q99_", type)) := vec_qevd(0.99, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV")
  )
}

# Overall by site
yearly_eva_params <- yearly_max_by_site %>%
  group_by(site, simulation_run) %>%
  summarise(
    across(all_of(starts_with("yearly_max_")), fit_fevd_and_get_params),
    .groups = "drop"
  )

yearly_eva_params <- tibble(do.call(data.frame, yearly_eva_params))

yearly_eva_params %<>% 
  get_quantiles_from_params_of_data_type("yearly_max_real") %>%
  get_quantiles_from_params_of_data_type("yearly_max_galm") %>%
  get_quantiles_from_params_of_data_type("yearly_max_bamlss")


yearly_eva_params %>% dplyr::select(
  site, 
  simulation_run,
  yearly_max_real.shape,
  yearly_max_galm.shape,
  yearly_max_bamlss.shape,
  q50_yearly_max_real,
  q50_yearly_max_galm,
  q50_yearly_max_bamlss,
  q90_yearly_max_real,
  q90_yearly_max_galm,
  q90_yearly_max_bamlss,
  q99_yearly_max_real,
  q99_yearly_max_galm,
  q99_yearly_max_bamlss
)

yearly_eva_params %>% dplyr::select(
  site, 
  yearly_max_real.shape,
  yearly_max_galm.shape,
  yearly_max_bamlss.shape,
  simulation_run
) %>% group_by(site) %>% summarise(across(ends_with("shape"), ~mean(.x)), .groups = "drop") %>% filter(site == "pooled") 

# Mean average distance: shape
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(yearly_max_real.shape - yearly_max_galm.shape)),
  mean_average_dist_bamlss = mean(abs(yearly_max_real.shape - yearly_max_bamlss.shape))
)

# location
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(yearly_max_real.location - yearly_max_galm.location)),
  mean_average_dist_bamlss = mean(abs(yearly_max_real.location - yearly_max_bamlss.location))
)

# scale
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(yearly_max_real.scale - yearly_max_galm.scale)),
  mean_average_dist_bamlss = mean(abs(yearly_max_real.scale - yearly_max_bamlss.scale))
)

# q50
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(q50_yearly_max_real - q50_yearly_max_galm)),
  mean_average_dist_bamlss = mean(abs(q50_yearly_max_real - q50_yearly_max_bamlss))
)

# q90
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(q90_yearly_max_real - q90_yearly_max_galm)),
  mean_average_dist_bamlss = mean(abs(q90_yearly_max_real - q90_yearly_max_bamlss))
)

# q99
yearly_eva_params %>% summarise(
  mean_average_dist_galm = mean(abs(q99_yearly_max_real - q99_yearly_max_galm)),
  mean_average_dist_bamlss = mean(abs(q99_yearly_max_real - q99_yearly_max_bamlss))
)

# GEV Analysis -- Seasonal maxima --------------------------------------------------------------


seasonal_max_by_site <- results %>%
  filter(site %in% sites_for_gev) %>%
  group_by(
    year,
    meteorological_season,
    site,
    simulation_run
  ) %>%
  summarise(
    seasonal_max_real = max(real_precip),
    seasonal_max_galm = max(simulation_values_galm),
    seasonal_max_bamlss = max(simulation_values_bamlss),
    .groups = "drop"
  )

seasonal_max_by_site <- rbind(seasonal_max_by_site, seasonal_max_by_site %>% mutate(site = "pooled"))

library("extRemes")

vec_qevd <- Vectorize(qevd)

fit_fevd_and_get_params <- function(data) {
  fit <- fevd(data, type = "GEV", method = "MLE")
  return(tibble(
    location = distill.fevd(fit)["location"],
    scale = distill.fevd(fit)["scale"],
    shape = distill.fevd(fit)["shape"]
  ))
}

get_quantiles_from_params_of_data_type <- function(.data, type) {
  mutate(
    .data,
    !!as.symbol(paste0("q50_", type)) := vec_qevd(0.5, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV"),
    !!as.symbol(paste0("q90_", type)) := vec_qevd(0.9, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV"),
    !!as.symbol(paste0("q99_", type)) := vec_qevd(0.99, loc = !!as.symbol(paste0(type, ".location")), scale = !!as.symbol(paste0(type, ".scale")), shape = !!as.symbol(paste0(type, ".shape")), type = "GEV")
  )
}

seasonal_eva_params <- seasonal_max_by_site %>%
  group_by(site, meteorological_season, simulation_run) %>%
  summarise(
    across(all_of(starts_with("seasonal_max_")), fit_fevd_and_get_params),
    .groups = "drop"
  )

seasonal_eva_params <- tibble(do.call(data.frame, seasonal_eva_params))

seasonal_eva_params %<>% 
  get_quantiles_from_params_of_data_type("seasonal_max_real") %>%
  get_quantiles_from_params_of_data_type("seasonal_max_galm") %>%
  get_quantiles_from_params_of_data_type("seasonal_max_bamlss")


seasonal_eva_params %>% dplyr::select(
  site, 
  simulation_run,
  meteorological_season,
  seasonal_max_real.shape,
  seasonal_max_galm.shape,
  seasonal_max_bamlss.shape,
  q50_seasonal_max_real,
  q50_seasonal_max_galm,
  q50_seasonal_max_bamlss,
  q90_seasonal_max_real,
  q90_seasonal_max_galm,
  q90_seasonal_max_bamlss,
  q99_seasonal_max_real,
  q99_seasonal_max_galm,
  q99_seasonal_max_bamlss
)

seasonal_eva_params %>% dplyr::select(
  site, 
  simulation_run,
  meteorological_season,
  seasonal_max_real.shape,
  seasonal_max_galm.shape,
  seasonal_max_bamlss.shape,
) %>% group_by(site, meteorological_season) %>% summarise(across(ends_with("shape"), ~sd(.x)), .groups = "drop") %>% filter(site == "pooled") 


# Mean average distance: shape
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(seasonal_max_real.shape - seasonal_max_galm.shape)),
  mean_average_dist_bamlss = mean(abs(seasonal_max_real.shape - seasonal_max_bamlss.shape))
)

# location
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(seasonal_max_real.location - seasonal_max_galm.location)),
  mean_average_dist_bamlss = mean(abs(seasonal_max_real.location - seasonal_max_bamlss.location))
)

# scale
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(seasonal_max_real.scale - seasonal_max_galm.scale)),
  mean_average_dist_bamlss = mean(abs(seasonal_max_real.scale - seasonal_max_bamlss.scale))
)

# q50
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(q50_seasonal_max_real - q50_seasonal_max_galm)),
  mean_average_dist_bamlss = mean(abs(q50_seasonal_max_real - q50_seasonal_max_bamlss))
)

# q90
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(q90_seasonal_max_real - q90_seasonal_max_galm)),
  mean_average_dist_bamlss = mean(abs(q90_seasonal_max_real - q90_seasonal_max_bamlss))
)

# q99
seasonal_eva_params %>% group_by(meteorological_season) %>% summarise(
  mean_average_dist_galm = mean(abs(q99_seasonal_max_real - q99_seasonal_max_galm)),
  mean_average_dist_bamlss = mean(abs(q99_seasonal_max_real - q99_seasonal_max_bamlss))
)



# GEV Analysis -- Yearly maxima of average rainfall --------------------------------------------------------------

yearly_max_by_site <- results %>%
  group_by(date, year, simulation_run) %>%
  summarise(real_precip = mean(real_precip), simulation_values_galm = mean(simulation_values_galm), simulation_values_bamlss = mean(simulation_values_bamlss), .groups = "drop") %>%
  group_by(
    year,
    simulation_run
  ) %>%
  summarise(
    yearly_max_real = max(real_precip),
    yearly_max_galm = max(simulation_values_galm),
    yearly_max_bamlss = max(simulation_values_bamlss),
    .groups = "drop"
  )

# Overall by site
yearly_eva_params <- yearly_max_by_site %>%
  summarise(
    across(all_of(starts_with("yearly_max_")), fit_fevd_and_get_params),
    .groups = "drop"
  )

yearly_eva_params <- tibble(do.call(data.frame, yearly_eva_params))

yearly_eva_params %<>% 
  get_quantiles_from_params_of_data_type("yearly_max_real") %>%
  get_quantiles_from_params_of_data_type("yearly_max_galm") %>%
  get_quantiles_from_params_of_data_type("yearly_max_bamlss")


yearly_eva_params %>% dplyr::select(
  yearly_max_real.shape,
  yearly_max_galm.shape,
  yearly_max_bamlss.shape,
  q50_yearly_max_real,
  q50_yearly_max_galm,
  q50_yearly_max_bamlss,
  q90_yearly_max_real,
  q90_yearly_max_galm,
  q90_yearly_max_bamlss,
  q99_yearly_max_real,
  q99_yearly_max_galm,
  q99_yearly_max_bamlss
)


