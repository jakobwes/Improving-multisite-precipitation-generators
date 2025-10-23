
# 1. Computation in the gamlss package ------------------------------------

# Compute composite likelihood AIC for model evaluation
get_aic_composite_likelihood <- function(fit) {

  # vcov(fit, robust = TRUE) returns I^-1 K I^-1, which is the inverse of the Godambe information G = I K^-1 I (or G = H J^-1 H using the notation from Varin et al.)
  inverse_godambe_information <- vcov(fit, robust = TRUE)

  # vcov(fit, robust = FALSE) returns I^-1, which is the inverse of the observed information
  observed_information <- solve(vcov(fit, robust = FALSE))

  return(
    -2 * logLik(fit) + sum(diag(observed_information %*% inverse_godambe_information))
  )
}

print_aic_composite_likelihood <- function(fit) {
  aic <- as.numeric(get_aic_composite_likelihood(fit))

  return(
    paste0("AIC: ", prettyNum(round(aic, 1), big.mark = ","))
  )
}


# 2. Computation for glms using the sandwich package ----------------------

library("sandwich")

get_sandwich_aic_adj <- function(fit){
  M <- meat(fit)
  B <- bread(fit)
  
  return(
    -2 * logLik(fit) + sum(diag(M %*% B))
  )
}

print_sandwich_aic_adj <- function(fit) {
  aic <- as.numeric(get_sandwich_aic_adj(fit))
  
  return(
    paste0("AIC: ", prettyNum(round(aic, 1), big.mark = ","))
  )
}
