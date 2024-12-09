library(rstan)
#' Title
#' stan_fit
#'
#' @param model_path         # a path string to the stan file
#' @param stan_data         # a list of data aligning to the formulation in the stan code used to take synthetic data into rstan
#'
#'
#' @return Explain return
#' A stan fit object including the sampling results is returned
#' @export
#'
#' @examples

stan_fit <- function(model_path = "/hmm.stan", stan_data, num_chains = 4, iter_times = 2000, warm_up = 1000){
  fit <- sampling(stan_model(model_code = readLines(model_path)),
                  data = stan_data,
                  chains = num_chains,
                  iter = iter_times,
                  warmup = warm_up,
                  #init = rep(list(init_fun()), 4),
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  return(fit)
}




plot_densities <- function(fit, params) {
  # Extract samples
  samples <- rstan::extract(fit, params)

  # Calculate number of rows and columns for plot layout
  n_params <- length(params)
  n_cols <- min(3, n_params)  # maximum 3 columns
  n_rows <- ceiling(n_params/n_cols)

  # Set up plot layout with reasonable margins
  old_par <- par(no.readonly = TRUE)  # save old parameters
  par(mfrow = c(n_rows, n_cols),
      mar = c(4, 4, 2, 1),  # reduce margins
      oma = c(0, 0, 2, 0))  # reduce outer margins

  # Create density plots
  for(param in params) {
    density_obj <- density(samples[[param]])
    plot(density_obj, main=param, xlab="Value", ylab="Density")
  }

  # Reset plotting parameters
  par(old_par)
}




# Posterior Predictive Checks
posterior_predictive_check <- function(fit, y_obs, n_samples = 100) {
  # Extract model parameters
  A_samples <- rstan::extract(fit, "A")$A      # transition probabilities
  mu_samples <- rstan::extract(fit, "mu")$mu   # means
  sigma_samples <- rstan::extract(fit, "sigma")$sigma  # standard deviations

  # Generate simulated data
  n_obs <- length(y_obs)
  y_rep <- matrix(0, n_samples, n_obs)

  for(i in 1:n_samples) {
    # Sample random iteration from posterior
    iter <- sample(1:nrow(A_samples), 1)

    # Generate new sequence using parameters from this iteration
    states <- numeric(n_obs)
    y_rep[i,] <- numeric(n_obs)

    # Initial state
    states[1] <- sample(1:ncol(A_samples[iter,,]), 1)
    y_rep[i,1] <- rnorm(1, mu_samples[iter, states[1]], sigma_samples[iter, states[1]])

    # Generate rest of sequence
    for(t in 2:n_obs) {
      states[t] <- sample(1:ncol(A_samples[iter,,]), 1,
                          prob = A_samples[iter, states[t-1],])
      y_rep[i,t] <- rnorm(1, mu_samples[iter, states[t]],
                          sigma_samples[iter, states[t]])
    }
  }

  # Plot
  plot(density(y_obs), main="Posterior Predictive Check",
       xlab="Value", ylab="Density", lwd=2)

  # Add simulated densities
  for(i in 1:n_samples) {
    lines(density(y_rep[i,]), col=rgb(0,0,1,0.1))
  }
  lines(density(y_obs), lwd=2)
  legend("topright", c("Observed", "Simulated"),
         lty=1, col=c("black", rgb(0,0,1,0.5)))
}

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



# Example usage:
# fit <- your stan fit object
# y_obs <- your observed data
# params <- c("mu", "sigma", "theta")  # your parameters of interest
# evaluate_fit(fit, y_obs, params)


