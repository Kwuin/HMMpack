# Function to generate initial values
source("R/sampler.R")
source("R/models.R")
init_fun <- function() {
  list(
    pi = rep(1/K, K),
    A = matrix(1/K, K, K),
    mu = sort(rnorm(K, mean=mean(stan_data$y), sd=sd(stan_data$y))),
    sigma = rep(sd(stan_data$y), K)
  )
}

# Set seed for reproducibility
seed = 123
set.seed(seed)

# Simulate data
N <- 200  # shorter sequence for testing
K <- 3    # number of states

# True parameters (more separated means for better identification)
pi_true <- c(0.6, 0.3, 0.1)
A_true <- matrix(c(0.8, 0.1, 0.1,
                   0.1, 0.8, 0.1,
                   0.1, 0.1, 0.8),
                 nrow = K, byrow = TRUE)
mu_true <- c(-3, 0, 3)
sigma_true <- c(0.5, 0.5, 0.5)

hmm_model = HMM(K,A_true, pi_true, seed = seed)
hmm_gaussian_model = HMM_Gaussian_Model(hmm_model, sigma_true^2)


y <- hmm_gaussian_model$generate_gaussian_observations(N)$observations

stan_data <- list(
  N = N,
  K = K,
  y = y
)

fit = stan_fit("hmm.stan", stan_data)


library(rstan)
library(bayesplot)
library(ggplot2)


matrix_param_names <- c(outer(1:K, 1:K, FUN=function(i,j) paste0("A[", i, ",", j, "]")))
params = matrix_param_names

posterior_samples <- rstan::extract(fit, matrix_param_names)

summary <- summary(fit, matrix_param_names)$summary
# Add more informative column names
colnames(summary) <- c("Mean", "SE", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "N_eff", "Rhat")

rhat <- summary(fit)$summary[, "Rhat"]
problems <- which(rhat > 1.1)

# Check effective sample size
neff <- summary(fit)$summary[, "n_eff"]
low_neff <- which(neff < 100)

rstan::traceplot(fit, params)

neff_ratio(fit, params)

mcmc_acf(fit, params)

source("R/sampler.R")
plot_densities(fit, params)

posterior_predictive_check(fit, y, N)




