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


# 4. Define parameters to evaluate
# Create names for transition matrix parameters
matrix_param_names <- c(outer(1:K, 1:K, FUN=function(i,j) paste0("A[", i, ",", j, "]")))
# Add other parameters
param_names <- c(matrix_param_names,
                 paste0("mu[", 1:K, "]"),
                 paste0("sigma[", 1:K, "]"))

# 5. Run evaluation
results <- evaluate_stan_fit(fit, param_names, y, N, K)

# 6. Check results
# Look at summary statistics
print("Parameter Summaries:")
print(results$summary)

# Check for convergence issues
if(length(results$rhat_problems) > 0) {
  print("Parameters with convergence issues:")
  print(param_names[results$rhat_problems])
}

if(length(results$neff_problems) > 0) {
  print("Parameters with low effective sample size:")
  print(param_names[results$neff_problems])
}

# Compare with true values
print("Comparing estimates with true values:")
# For transition matrix
for(i in 1:K) {
  for(j in 1:K) {
    true_val <- A_true[i,j]
    est_val <- results$summary[paste0("A[",i,",",j,"]"), "Mean"]
    cat(sprintf("A[%d,%d]: True = %.3f, Estimated = %.3f\n",
                i, j, true_val, est_val))
  }
}

# For means and standard deviations
for(i in 1:K) {
  cat(sprintf("mu[%d]: True = %.3f, Estimated = %.3f\n",
              i, mu_true[i], results$summary[paste0("mu[",i,"]"), "Mean"]))
  cat(sprintf("sigma[%d]: True = %.3f, Estimated = %.3f\n",
              i, sigma_true[i], results$summary[paste0("sigma[",i,"]"), "Mean"]))
}
