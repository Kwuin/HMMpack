data {
  int<lower=1> N;              // length of observations
  int<lower=1> K;              // number of hidden states
  vector[N] y;                 // observations
}

transformed data {
  vector[K] alpha = rep_vector(2.0, K);  // prior concentration
}

parameters {
  simplex[K] pi;               // initial state probabilities
  simplex[K] A[K];            // transition probabilities (K x K matrix)
  ordered[K] mu;               // means for emissions (ordered for identifiability)
  vector<lower=0.1>[K] sigma;  // standard deviations for emissions
}

transformed parameters {
  matrix[N, K] log_alpha;      // forward algorithm log probabilities

  // Forward algorithm
  for (k in 1:K)
    log_alpha[1, k] = log(pi[k]) + normal_lpdf(y[1] | mu[k], sigma[k]);

  for (t in 2:N) {
    for (j in 1:K) {
      vector[K] temp;
      for (k in 1:K)
        temp[k] = log_alpha[t-1, k] + log(A[k, j]);
      log_alpha[t, j] = log_sum_exp(temp) + normal_lpdf(y[t] | mu[j], sigma[j]);
    }
  }
}

model {
  // Priors
  pi ~ dirichlet(alpha);
  for (k in 1:K)
    A[k] ~ dirichlet(alpha);

  mu ~ normal(0, 10);
  sigma ~ normal(0, 5);

  target += log_sum_exp(log_alpha[N]);
}

generated quantities {
  array[N] int<lower=1,upper=K> states;
  matrix[N, K] log_beta;

  { // local scope
    matrix[N, K] posterior;
    vector[K] work;

    // Initialize backward pass
    log_beta[N,] = rep_vector(0.0, K)';

    // Compute posterior for last time step
    posterior[N] = log_alpha[N];
    posterior[N] = posterior[N] - log_sum_exp(posterior[N]); // normalize

    // Sample last state
    for (k in 1:K)
      work[k] = posterior[N, k];

    states[N] = categorical_logit_rng(work);

    // Backward pass
    for (t in 1:(N-1)) {
      int t_curr = N - t;
      int t_next = t_curr + 1;

      // Compute backward probabilities
      for (k in 1:K) {
        log_beta[t_curr, k] = log(A[k, states[t_next]]) +
                             normal_lpdf(y[t_next] | mu[states[t_next]], sigma[states[t_next]]) +
                             log_beta[t_next, states[t_next]];
      }

      // Compute posterior
      posterior[t_curr] = log_alpha[t_curr] + log_beta[t_curr];
      posterior[t_curr] = posterior[t_curr] - log_sum_exp(posterior[t_curr]); // normalize

      // Sample state
      for (k in 1:K)
        work[k] = posterior[t_curr, k];

      states[t_curr] = categorical_logit_rng(work);
    }
  }
}
