library("tidyverse")
library("magrittr")
library("lubridate")

# Structure
# 1. Global helper functions
# 3. Stratified bootstrap with days as strata
#     3.1 Bootstrap returning ids in rainy_days
#     3.2 Bootstrap returning weighted observation matrix for weighted likelihood analysis
# 4. Stratified bootstrap of sites
# 5. Double stratified sampling


# 1. Global helper functions ----------------------------------------------

# Helper function for stratified bootstrap:
get_idxs_for_bootstrap_within_one_group <- function(group_name, vector_of_groups, n = NULL) {
  # vector_of_groups: a list of group ids where the idxs correspond to the idxs of observations
  # and the group id named to the group this observation is in.

  ixds_in_group <- which(vector_of_groups == group_name)
  if (is.null(n)) n <- length(ixds_in_group)

  # Special case where strata only has one observation
  if (length(ixds_in_group) == 1) {
    return(rep(ixds_in_group, n))
  }

  return(sample(ixds_in_group, n, replace = TRUE))
}
get_idxs_for_stratified_bootstrap <- function(vector_of_groups) {
  idxs <- c()
  for (group in unique(vector_of_groups)) {
    idxs <- c(idxs, get_idxs_for_bootstrap_within_one_group(group, vector_of_groups))
  }
  return(idxs)
}


# 2. Stratified bootstrap -------------------------------------------------

# This is a bootstrapping strategy preserving entire days.
# Method: Sample days instead of individual observations with the goal of having same number of days and observations as in the initial sample
#
# For that:
#   1. Create a pool of all days in the given sample
#   2. Group together all days with the same number of non-missing observations as strata
#   3. Use stratified sampling to sample from this dataset, so sample from each strata as many items as the strata has in the initial sample, with replacement
#   4. Put the observations of the days in the final dataset

get_observations_corresponding_to_days_in_bootstrapped_dataset <- function(day_nr_of_days_in_bootstrap, rainy_days) {
  idxs_observations <- c()
  for (day_nr in day_nr_of_days_in_bootstrap) {
    idxs_observations <- c(idxs_observations, which(rainy_days$day_nr == day_nr))
  }
  return(idxs_observations)
}

# 2.1 Bootstrap returning ids in rainy_days

# rainy_days: usual dataframe of observations with day_nr (day index)
get_ids_for_stratified_bootstrap <- function(rainy_days, n_bootstraps = 25) {
  bootstrapped_dataset_ids <- matrix(0, dim(rainy_days)[1], n_bootstraps)
  bootstrapped_day_nrs <- matrix(0, length(unique(rainy_days$day_nr)), n_bootstraps)

  # 1. Create a pool of all days in the given sample and 2. Group together all days with the same number of non-missing observations as strata
  temp_observations_per_day <- rainy_days %>% count(day_nr)

  # Add a unique identifier to each grouped together number of days
  temp_strata_nr_and_match_to_nr_of_observations <- data.frame(strata_number = c(1:length(unique(temp_observations_per_day$n))), observations_per_day_in_this_strata = unique(temp_observations_per_day$n))
  temp_day_nr_with_strata_nr <- temp_observations_per_day %>% left_join(temp_strata_nr_and_match_to_nr_of_observations, by = c("n" = "observations_per_day_in_this_strata"))

  for (temp_i in 1:n_bootstraps) {
    # 3. Use stratified sampling to sample from this dataset, so sample from each strata as many items as the strata has in the initial sample, with replacement
    temp_day_nr_of_days_in_bootstrap <- temp_day_nr_with_strata_nr$day_nr[get_idxs_for_stratified_bootstrap(temp_day_nr_with_strata_nr$strata_number)]
    # 4. Put the observations of the days in the final dataset
    bootstrapped_dataset_ids[, temp_i] <- get_observations_corresponding_to_days_in_bootstrapped_dataset(temp_day_nr_of_days_in_bootstrap, rainy_days)
    # Save the days resampled (for comparison)
    bootstrapped_day_nrs[, temp_i] <- temp_day_nr_of_days_in_bootstrap
  }
  return(list(
    bootstrapped_dataset_ids = bootstrapped_dataset_ids,
    bootstrapped_day_nrs = bootstrapped_day_nrs
  ))
}

# 2.2 Bootstrap returning weighted observation matrix for weighted likelihood analysis

get_weights_matrix_for_stratified_bootstrap <- function(rainy_days, n_bootstraps = 25) {
  temp_stratified_sample <- get_ids_for_stratified_bootstrap(rainy_days, n_bootstraps)

  bootstrapped_dataset_ids <- temp_stratified_sample$bootstrapped_dataset_ids
  bootstrapped_day_nrs <- temp_stratified_sample$bootstrapped_day_nrs

  bootstrapped_weights_matrix <- matrix(0, dim(bootstrapped_dataset_ids)[1], dim(bootstrapped_dataset_ids)[2])

  for (i in 1:dim(bootstrapped_weights_matrix)[1]) {
    for (j in 1:dim(bootstrapped_weights_matrix)[2]) {
      bootstrapped_weights_matrix[bootstrapped_dataset_ids[i, j], j] <- bootstrapped_weights_matrix[bootstrapped_dataset_ids[i, j], j] + 1
    }
  }

  return(list(
    bootstrapped_weights_matrix = bootstrapped_weights_matrix,
    bootstrapped_day_nrs = bootstrapped_day_nrs
  ))
}


# 3. Stratified bootstrap of sites -------------------------------------------

get_ids_for_stratified_bootstrap_of_sites <- function(rainy_days) {
  return(get_idxs_for_stratified_bootstrap(rainy_days$site))
}

get_weights_matrix_for_stratified_bootstrap_of_sites <- function(groups, B) {
  n <- length(groups)
  folds <- matrix(0, nrow = n, ncol = B)

  for (s in levels(groups)) {
    indx <- which(groups == s)
    n_i <- length(indx)

    folds[indx, ] <- rmultinom(B, n_i, rep(1, n_i))
  }
  return(folds)
}


# 4. Double stratified sampling -------------------------------------------

# This is a bootstrapping strategy preserving entire days, whilst also preserving the amount of rainy sites on each day.
# Method: Group together days with same number of observations recorded and same number of rainy sites. Do stratified bootstrap of days, with these strata.

# This is a double stratified bootstrap, since the marginal (when only looking at the rainy sites) is a stratified bootstrap, but also the full (when only looking at non missing observations)

get_ids_for_double_stratified_bootstrap <- function(dat, n_bootstraps = 25) {
  # dat: usual dataframe of observations with day_nr (day index)

  bootstrapped_dataset_ids <- matrix(0, dim(dat)[1], n_bootstraps)
  bootstrapped_day_nrs <- matrix(0, length(unique(dat$day_nr)), n_bootstraps)

  rainy_days <- dat %>% filter(precip_occurrence == 1)

  bootstrapped_dataset_ids_in_rainy_days <- matrix(0, dim(rainy_days)[1], n_bootstraps)

  # 1. Create a pool of all days in the given sample and 2. Group together all days with the same number of non-missing observations and rainy sites as strata
  day_nrs_with_grouping_variables <- dat %>%
    group_by(day_nr) %>%
    summarise(n_observations = n(), n_rainy_sites = sum(precip_occurrence))

  # Defines the strata (groups) and gives them a group-id
  strata_and_defining_categories <- day_nrs_with_grouping_variables %>%
    group_by(n_observations, n_rainy_sites) %>%
    summarise(n_days_in_strata = n(), .groups = "drop")
  strata_and_defining_categories$strata_id <- 1:dim(strata_and_defining_categories)[1]

  day_nrs_with_grouping_variables <- day_nrs_with_grouping_variables %>% left_join(strata_and_defining_categories, by = c("n_observations" = "n_observations", "n_rainy_sites" = "n_rainy_sites"))

  for (temp_i in 1:n_bootstraps) {
    # 3. Use stratified sampling to sample from this dataset, so sample from each strata as many items as the strata has in the initial sample, with replacement
    day_nr_of_days_in_bootstrap <- day_nrs_with_grouping_variables$day_nr[get_idxs_for_stratified_bootstrap(factor(day_nrs_with_grouping_variables$strata_id))]
    # 4. Put the observations of the days in the final dataset
    bootstrapped_dataset_ids[, temp_i] <- get_observations_corresponding_to_days_in_bootstrapped_dataset(day_nr_of_days_in_bootstrap, dat)
    # Save the days resampled (for comparison)
    bootstrapped_day_nrs[, temp_i] <- day_nr_of_days_in_bootstrap
    # Get idxs in rainy_days for amounts model
    bootstrapped_dataset_ids_in_rainy_days[, temp_i] <- get_observations_corresponding_to_days_in_bootstrapped_dataset(day_nr_of_days_in_bootstrap, rainy_days)
  }
  return(list(
    bootstrapped_dataset_ids = bootstrapped_dataset_ids,
    bootstrapped_day_nrs = bootstrapped_day_nrs,
    bootstrapped_dataset_ids_in_rainy_days = bootstrapped_dataset_ids_in_rainy_days
  ))
}



# 6. Apply ----------------------------------------------------------------

# dat <- read_csv("data/3_data_processed/processed_data_covariates_only.csv") %>% 
#   mutate(
#     precip_occurrence = if_else(precip > 0, 1, 0),
#     day_nr = as.integer(date - min(date))
#   )
# 
# # Get double stratified bootstrap
# double_stratified_bootstrap <- get_ids_for_double_stratified_bootstrap(dat)
# 
# write.table(double_stratified_bootstrap$bootstrapped_dataset_ids, file="src/2_model_selection/2_nonparametric_model_selection/0_bootstrap/bootstrapped_dataset_ids.txt", row.names=FALSE, col.names=FALSE)
# write.table(double_stratified_bootstrap$bootstrapped_day_nrs, file="src/2_model_selection/2_nonparametric_model_selection/0_bootstrap/bootstrapped_day_nrs.txt", row.names=FALSE, col.names=FALSE)
# write.table(double_stratified_bootstrap$bootstrapped_dataset_ids_in_rainy_days, file="src/2_model_selection/2_nonparametric_model_selection/0_bootstrap/bootstrapped_dataset_ids_in_rainy_days.txt", row.names=FALSE, col.names=FALSE)
