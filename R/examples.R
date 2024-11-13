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
mu_true <- c(-3, 0, 3)
sigma_true <- c(0.5, 0.5, 0.5)

hmm_model = HMM(K, pi_true, A_true)
hmm_gaussian_model = HMM_Gaussian_Model(hmm_model, sigma_true^2)
