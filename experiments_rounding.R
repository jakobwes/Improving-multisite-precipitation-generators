rm(list = ls())
library("gamlss")
library("tidyverse")
library("lubridate")
library("gamlss")

make_qq_plot <- function(y, sim_y){
  stats::qqplot(
    y, sim_y,
    xlab = "Observed Values",
    ylab = "Simulated Values",
    pch = 20,      # solid points
    cex = 1.0      # smaller size
  )
  
  # add 1–1 line
  abline(0, 1, col = "red", lty = 2, lwd = 2)
}

# 1. Simulation of rounding only  -----------------------------------------------------------

# Unrounded values
x <- rGA(100000, 5, 1)

# Rounded values
y <- round(x*2)/2
y <- y[y > 0]

# Fit model
fit <- gamlss(formula = y ~ 1, family = GA(mu.link = "log", sigma.link = "log"))

# Plot density of unrounded vs. rounded values
plot(density(x), col = "red"); lines(density(y), col = "blue")

# Simulate values
sim_y <- rGA(length(y), predict(fit, what = "mu", type = "response"),  predict(fit, what = "sigma", type = "response"))

# Plot density of rounded vs. simulated values
plot(density(y), col = "red"); lines(density(sim_y), col = "blue")

png("plots/rounded_vs_simulated.png", width = 400, height = 350)
make_qq_plot(y, sim_y)
dev.off()


# Plot density of rounded vs. rounded simulated values
sim_y <- round(sim_y*2)/2
sim_y <- sim_y[sim_y > 0]
plot(density(y), col = "red"); lines(density(sim_y), col = "blue")
png("plots/rounded_vs_rounded_simulated.png", width = 400, height = 350)
make_qq_plot(y, sim_y)
dev.off()


# 2. Simulation of soft thresholding  -----------------------------------------------------------

# Unrounded values
x <- rGA(100000, 5, 1)

# Rounded values
y <- round(x*2)/2
y <- y[y > 0]
y <- y - 0.49

# Fit model
fit <- gamlss(formula = y ~ 1, family = GA(mu.link = "log", sigma.link = "log"))

# Plot density of unrounded vs. rounded values
plot(density(x), col = "red"); lines(density(y), col = "blue")

# Simulate values
sim_y <- rGA(length(y), predict(fit, what = "mu", type = "response"),  predict(fit, what = "sigma", type = "response"))

# Plot density of rounded vs. simulated values
plot(density(y), col = "red"); lines(density(sim_y), col = "blue")

png("plots/rounded_vs_simulated_soft_thresholding.png", width = 400, height = 350)
make_qq_plot(y, sim_y)
dev.off()


# Plot density of rounded vs. rounded simulated values
sim_y <- round(sim_y*2)/2
sim_y <- sim_y[sim_y > 0]
sim_y <- sim_y - 0.49

plot(density(y), col = "red"); lines(density(sim_y), col = "blue")
png("plots/rounded_vs_rounded_simulated_soft_thresholding.png", width = 400, height = 350)
make_qq_plot(y, sim_y)
dev.off()


# CRPS-based optimisation
library('scoringRules')

optim_gamma_crps <- function(x){
  return(mean(crps_gamma(y, shape = exp(x[1]), rate = exp(x[2]))))
}
optim_res <- optim(c(0,0), optim_gamma_crps, method = 'BFGS')
sim_y <- rgamma(length(y), shape = exp(optim_res$par[1]), rate = exp(optim_res$par[2]))
make_qq_plot(y, sim_y)


# QQplot of residuals
sim_y <- round(x*10)/10
sim_y <- sim_y[sim_y > 0]
sim_y <- sim_y - 0.09
fit <- gamlss(formula = y ~ 1, family = GA(mu.link = "log", sigma.link = "log"))
qq <- qnorm(pGA(y, predict(fit, what = "mu", type = "response"),  predict(fit, what = "sigma", type = "response")))
qqnorm(qq); qqline(qq)

