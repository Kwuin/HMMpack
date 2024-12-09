# Hidden Markov Model for sequential data
This is a R package for implementing Hidden Markov Models(HMM) for sequential data. It includes functionalities of generating synthetic data with arbitray transition matrices, initial distribution and differnet emission distribution for continuous data. It can also be used to do Bayesian inference using Rstan to do Hamiltonian Monte Carlo sampling method. Calculation of log likelihood for given sequential data and model is provided. 

## Installation
Run the following code in Rstudio: 
```{r}
devtools::install_github("Kwuin/HMMpack")
```
## small example 
```{r}
library(HMMpack)

# Generate sample data
N <- 100  # sequence length
K <- 3    # number of states
A <- matrix(c(0.8, 0.1, 0.1,
              0.1, 0.8, 0.1,
              0.1, 0.1, 0.8), nrow=K)
mu <- c(-2, 0, 2)
sigma <- c(0.5, 0.5, 0.5)

# Simulate data
sim_data <- simulate_hmm(N, K, A, mu, sigma)
```
# Fit model
```{r}
fit <- fit_hmm(sim_data, K)
```


