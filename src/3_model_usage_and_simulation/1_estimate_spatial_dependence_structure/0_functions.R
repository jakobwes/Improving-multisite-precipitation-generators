# Functions
library(mvtnorm)

get_pivoted_day_classes_for_sites <- function(dat_with_a_b, site1, site2) {
  
  # Choose dates on which both sites haves measurements
  dat_specific_sites <- dat_with_a_b %>%
    filter(site == site1 | site == site2) %>%
    group_by(day_nr) %>%
    filter(dplyr::n() == 2) %>%
    ungroup()

  if (nrow(dat_specific_sites) == 0) {
    days_both_wet <- data.frame()
    days_site1_wet_site2_dry <- data.frame()
    days_site2_wet_site1_dry <- data.frame()
    days_both_dry <- data.frame()
  } else {
    dat_specific_sites %<>%
      dplyr::select(
        site,
        date,
        precip,
        precip_occurrence,
        a,
        b
      ) %>%
      pivot_wider(names_from = site, values_from = c(precip, precip_occurrence, a, b))

    days_both_wet <- dat_specific_sites %>% filter(!!as.symbol(paste0("precip_occurrence_", site1)) == TRUE & !!as.symbol(paste0("precip_occurrence_", site2)) == TRUE)
    days_site1_wet_site2_dry <- dat_specific_sites %>% filter(!!as.symbol(paste0("precip_occurrence_", site1)) == TRUE & !!as.symbol(paste0("precip_occurrence_", site2)) == FALSE)
    days_site2_wet_site1_dry <- dat_specific_sites %>% filter(!!as.symbol(paste0("precip_occurrence_", site1)) == FALSE & !!as.symbol(paste0("precip_occurrence_", site2)) == TRUE)
    days_both_dry <- dat_specific_sites %>% filter(!!as.symbol(paste0("precip_occurrence_", site1)) == FALSE & !!as.symbol(paste0("precip_occurrence_", site2)) == FALSE)
  }

  return(
    list(
      days_both_wet = days_both_wet,
      days_site1_wet_site2_dry = days_site1_wet_site2_dry,
      days_site2_wet_site1_dry = days_site2_wet_site1_dry,
      days_both_dry = days_both_dry
    )
  )
}

get_log_likelihood <- function(rho, days_both_wet, days_site1_wet_site2_dry, days_site2_wet_site1_dry, days_both_dry) {
  cov_mat <- matrix(c(1, rho, rho, 1), 2, 2)

  # Site 1 wet, site 2 dry:
  if (nrow(days_site1_wet_site2_dry) == 0) {
    ll_site1_wet_site2_dry <- 0
  } else {
    b1 <- days_site1_wet_site2_dry %>% pull(!!as.symbol(paste0("b_", site1)))
    a2 <- days_site1_wet_site2_dry %>% pull(!!as.symbol(paste0("a_", site2)))

    ll_site1_wet_site2_dry <- sum(log(pnorm(a2, mean = rho * b1, sd = sqrt(1 - rho^2))))
  }

  # Site 2 wet, site 1 dry:
  if (nrow(days_site2_wet_site1_dry) == 0) {
    ll_site2_wet_site1_dry <- 0
  } else {
    b2 <- days_site2_wet_site1_dry %>% pull(!!as.symbol(paste0("b_", site2)))
    a1 <- days_site2_wet_site1_dry %>% pull(!!as.symbol(paste0("a_", site1)))

    ll_site2_wet_site1_dry <- sum(log(pnorm(a1, mean = rho * b2, sd = sqrt(1 - rho^2))))
  }

  # Both sites dry
  if (nrow(days_both_dry) == 0) {
    ll_both_dry <- 0
  } else {
    a1 <- days_both_dry %>% pull(!!as.symbol(paste0("a_", site1)))
    a2 <- days_both_dry %>% pull(!!as.symbol(paste0("a_", site2)))

    ll_both_dry <- sum(sapply(1:length(a1), function(i) log(pmvnorm(upper = c(a1[i], a2[i]), mean = c(0, 0), sigma = cov_mat))))
  }

  # Both sites wet
  if (nrow(days_both_wet) == 0) {
    ll_both_wet <- 0
  } else {
    b1 <- days_both_wet %>% pull(!!as.symbol(paste0("b_", site1)))
    b2 <- days_both_wet %>% pull(!!as.symbol(paste0("b_", site2)))

    ll_both_wet <- sum(sapply(1:length(b1), function(i) log(dmvnorm(x = c(b1[i], b2[i]), mean = c(0, 0), sigma = cov_mat))))
  }

  return(
    ll_site1_wet_site2_dry + ll_site2_wet_site1_dry + ll_both_dry + ll_both_wet
  )
}


# Matérn infilling

get_dist_mat <- function(sites, dat){
  
  nr_of_sites <- length(sites)
  dist_mat <- matrix(1, nr_of_sites, nr_of_sites)
  
  for (i in 1:nr_of_sites) {
    for (j in 1:i) {
      if (i == j) {
        next
      }
      
      site1 <- sites[i]
      site2 <- sites[j]
      
      lat_site_1 <- (dat %>% filter(site == site1) %>% pull(lat))[1]
      lon_site_1 <- (dat %>% filter(site == site1) %>% pull(lon))[1]
      
      lat_site_2 <- (dat %>% filter(site == site2) %>% pull(lat))[1]
      lon_site_2 <- (dat %>% filter(site == site2) %>% pull(lon))[1]
      
      dist <- sqrt((lat_site_1 - lat_site_2)^2 + (lon_site_1 - lon_site_2)^2)
      dist_mat[i,j] <- dist
      dist_mat[j,i] <- dist
    }
  }
  
  return(dist_mat)
}

get_matern <- function(rho, sigma, d, nu = 2){
  return(
    ((sigma^2 * 2^(1-nu))/gamma(nu)) * (sqrt(2*nu)*d/rho)^nu * besselK(2*d/rho, nu = 2)
  )
}

get_exponential <- function(rho, sigma, d, nu = 2){
  return(
    sigma^2 * exp(-d/rho)
  )
}

cost_function_matern <- function(rho, sigma, cov_mat, dist_mat, nu = 2){
  diag(cov_mat) <- NA
  covariances <- cov_mat[!is.na(cov_mat)]
  distances <- dist_mat[!is.na(cov_mat)]
  cost <- (get_matern(rho, sigma, distances, nu) - covariances)^2
  return(sum(cost))
}

cost_function_matern_weighted <- function(rho, sigma, cov_mat, dist_mat, n_datapoints_mat, nu = 2){
  diag(cov_mat) <- NA
  covariances <- cov_mat[!is.na(cov_mat)]
  distances <- dist_mat[!is.na(cov_mat)]
  N <- (sum(n_datapoints_mat)-nrow(n_datapoints_mat))/2 # minus as diag is filled with ones
  weights <- n_datapoints_mat[!is.na(cov_mat)]/N 
  cost <- (get_matern(rho, sigma, distances, nu) - covariances)^2
  return(sum(weights*cost))
}

cost_function_matern_weighted_all_params <- function(params, cov_mat, dist_mat, n_datapoints_mat){
  rho <- exp(params[1])
  sigma <- params[2]
  nu <- exp(params[3])
  
  diag(cov_mat) <- NA
  covariances <- cov_mat[!is.na(cov_mat)]
  distances <- dist_mat[!is.na(cov_mat)]
  N <- (sum(n_datapoints_mat)-nrow(n_datapoints_mat))/2 # minus as diag is filled with ones
  weights <- n_datapoints_mat[!is.na(cov_mat)]/N 
  cost <- (get_matern(rho, sigma, distances, nu) - covariances)^2
  return(sum(weights*cost))
}
