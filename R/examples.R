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


# Function to evaluate MCMC sampling results from Stan
evaluate_stan_fit <- function(fit, params, y, N, K) {
  # 1. Extract parameter summaries
  # Get posterior summary statistics for each parameter
  fit_summary <- rstan::summary(fit, pars = params)$summary
  colnames(fit_summary) <- c("Mean", "SE", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "N_eff", "Rhat")
  print("Parameter Summaries:")
  print(fit_summary)

  # 2. Convergence Diagnostics
  # Check R-hat values (should be close to 1)
  rhat <- fit_summary[, "Rhat"]
  rhat_problems <- which(rhat > 1.1)
  if(length(rhat_problems) > 0) {
    cat("\nParameters with high R-hat (>1.1):\n")
    print(params[rhat_problems])
  }

  # Check effective sample size (should be >100)
  neff <- fit_summary[, "N_eff"]
  low_neff_problems <- which(neff < 100)
  if(length(low_neff_problems) > 0) {
    cat("\nParameters with low effective sample size (<100):\n")
    print(params[low_neff_problems])
  }

  # 3. Visual Diagnostics
  # Set up plotting area
  par(mfrow=c(2,2))

  # Trace plots for checking mixing and convergence
  print("Generating trace plots...")
  rstan::traceplot(fit, pars = params)

  # Autocorrelation plots
  print("Generating autocorrelation plots...")
  if(require(bayesplot)) {
    mcmc_acf(fit, pars = params)
  }

  # Posterior density plots
  print("Generating density plots...")
  samples <- rstan::extract(fit, pars = params)
  for(param in params) {
    plot(density(samples[[param]]),
         main=paste("Density of", param),
         xlab=param)
  }

  # 4. Posterior Predictive Check
  print("Generating posterior predictive check...")
  # Extract parameters needed for simulation
  A_samples <- rstan::extract(fit, "A")$A
  mu_samples <- rstan::extract(fit, "mu")$mu
  sigma_samples <- rstan::extract(fit, "sigma")$sigma

  # Generate predictions
  n_pred <- 100  # number of predictions to generate
  y_pred <- matrix(0, n_pred, N)

  for(i in 1:n_pred) {
    # Sample one set of parameters
    idx <- sample(1:nrow(A_samples), 1)
    A <- A_samples[idx,,]
    mu <- mu_samples[idx,]
    sigma <- sigma_samples[idx,]

    # Generate sequence
    states <- numeric(N)
    states[1] <- sample(1:K, 1)
    y_pred[i,1] <- rnorm(1, mu[states[1]], sigma[states[1]])

    for(t in 2:N) {
      states[t] <- sample(1:K, 1, prob=A[states[t-1],])
      y_pred[i,t] <- rnorm(1, mu[states[t]], sigma[states[t]])
    }
  }

  # Plot observed vs predicted distributions
  plot(density(y), main="Posterior Predictive Check",
       xlab="Value", ylab="Density", lwd=2)
  for(i in 1:n_pred) {
    lines(density(y_pred[i,]), col=rgb(0,0,1,0.1))
  }
  lines(density(y), lwd=2)
  legend("topright", c("Observed", "Predicted"),
         lty=1, col=c("black", "blue"))

  # Reset plotting parameters
  par(mfrow=c(1,1))

  # 5. Return summary statistics
  return(list(
    summary = fit_summary,
    rhat_problems = rhat_problems,
    neff_problems = low_neff_problems
  ))
}


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

plot_all_diagnostics <- function(fit, params, y, N, K) {
  # 1. Trace plots (in a new window)
  dev.new(width=10, height=8)
  print("Generating trace plots...")
  rstan::traceplot(fit, params)

  # 2. Effective sample size ratio plot (in a new window)
  dev.new(width=10, height=6)
  print("Generating effective sample size ratio plot...")
  neff_ratio(fit, params)

  # 3. Autocorrelation plots (in a new window)
  dev.new(width=10, height=8)
  print("Generating autocorrelation plots...")
  mcmc_acf(fit, params)

  # 4. Density plots (in a new window)
  dev.new(width=10, height=8)
  print("Generating density plots...")
  # Calculate number of rows and columns for density plots
  n_params <- length(params)
  n_cols <- min(3, n_params)
  n_rows <- ceiling(n_params/n_cols)
  par(mfrow=c(n_rows, n_cols), mar=c(4,4,2,1))
  plot_densities(fit, params)

  # 5. Posterior predictive check (in a new window)
  dev.new(width=8, height=6)
  print("Generating posterior predictive check...")
  posterior_predictive_check(fit, y, N)
}

