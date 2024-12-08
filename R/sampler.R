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



# Example usage:
# fit <- your stan fit object
# y_obs <- your observed data
# params <- c("mu", "sigma", "theta")  # your parameters of interest
# evaluate_fit(fit, y_obs, params)


