library(rstan)

# Function to generate initial values
init_fun <- function() {
  list(
    pi = rep(1/K, K),
    A = matrix(1/K, K, K),
    mu = sort(rnorm(K, mean=mean(stan_data$y), sd=sd(stan_data$y))),
    sigma = rep(sd(stan_data$y), K)
  )
}

# Set seed for reproducibility
set.seed(123)

# Simulate data
N <- 200  # shorter sequence for testing
K <- 3    # number of states

# True parameters (more separated means for better identification)
pi_true <- c(0.6, 0.3, 0.1)
A_true <- matrix(c(0.8, 0.1, 0.1,
                   0.1, 0.8, 0.1,
                   0.1, 0.1, 0.8),
                 nrow = K, byrow = TRUE)
mu_true <- c(-3, 0, 3)  # more separated means
sigma_true <- c(0.5, 0.5, 0.5)

# Simulate data
simulate_hmm_data <- function(N, K, pi, A, mu, sigma) {
  states <- numeric(N)
  observations <- numeric(N)

  states[1] <- sample(1:K, 1, prob = pi)
  observations[1] <- rnorm(1, mu[states[1]], sigma[states[1]])

  for(t in 2:N) {
    states[t] <- sample(1:K, 1, prob = A[states[t-1], ])
    observations[t] <- rnorm(1, mu[states[t]], sigma[states[t]])
  }

  list(states = states, observations = observations)
}

sim_data <- simulate_hmm_data(N, K, pi_true, A_true, mu_true, sigma_true)

# Prepare data for Stan
stan_data <- list(
  N = N,
  K = K,
  y = sim_data$observations
)

# Compile and fit with conservative settings
hmm_model <- stan_model(model_code = readLines("hmm.stan"))
fit <- sampling(hmm_model,
                data = stan_data,
                chains = 4,
                iter = 2000,
                warmup = 1000,
                init = rep(list(init_fun()), 4),  # Use our initialization function
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15))

# Check results
print(fit, pars = c("mu", "sigma", "pi"))

# Plot results
plot(sim_data$observations, type = "l", col = "gray")
points(sim_data$states, col = "red", pch = ".")
states_est <- extract(fit)$states[1,]
points(states_est, col = "blue", pch = ".")
