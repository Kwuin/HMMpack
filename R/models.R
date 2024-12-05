#' Title
#' HMM Class Definition
#'
#' @param n_states    # a integer stating the number of states of the latent variable, represented by K
#' @param trans_mat   # K * K transition matrix of the Markov Chain formed by latent states
#' @param init_prob   # a numeric vector of K, the starting probability of each states
#' @param state_names # a vector of K strings, the names of each state
#'
#' @return
#' This is a class of HMM models to store the features of a Hidden Markov Model including the starting probabilities,
#' the transition matrices and number of states. It also provides a method to generate sequences of latent variables
#' with given parameters. It provides a method to generate synthetic continuous data and also a log likelihood calculator
#' for given data, which can be useful if an EM algorithm or variational inference method is to be implemented.
#' @export
#'
#' @examples

HMM <- setRefClass("HMM",
                   fields = list(
                     n_states = "numeric",          # Number of hidden states
                     trans_mat = "matrix",          # Transition probability matrix
                     init_prob = "numeric",         # Initial state probabilities
                     state_names = "character",     # Names of the states (optional)
                     seed = "numeric"               #seed
                   ),


                   methods = list(
                     # Initialize the HMM
                     initialize = function(n_states = 10, trans_mat = NULL, init_prob = NULL, state_names = NULL, seed = 123) {
                       "Initialize HMM with number of states, transition matrix, and initial probabilities"
                       n_states <<- n_states

                       # Set default transition matrix if not provided
                       if (is.null(trans_mat)) {
                         trans_mat <<- matrix(1/n_states, nrow = n_states, ncol = n_states)
                       } else {
                         if (!all(dim(trans_mat) == c(n_states, n_states))) {
                           stop("Transition matrix dimensions must match number of states")
                         }
                         if (!all(abs(rowSums(trans_mat) - 1) < 1e-10)) {
                           stop("Transition matrix rows must sum to 1")
                         }
                         trans_mat <<- trans_mat
                       }

                       # Set default initial probabilities if not provided
                       if (is.null(init_prob)) {
                         init_prob <<- rep(1/n_states, n_states)
                       } else {
                         if (length(init_prob) != n_states) {
                           stop("Initial probability vector length must match number of states")
                         }
                         if (abs(sum(init_prob) - 1) > 1e-10) {
                           stop("Initial probabilities must sum to 1")
                         }
                         init_prob <<- init_prob
                       }

                       # Set state names if provided
                       if (!is.null(state_names)) {
                         if (length(state_names) != n_states) {
                           stop("Number of state names must match number of states")
                         }
                         state_names <<- state_names
                       } else {
                         state_names <<- paste0("State", 1:n_states)
                       }

                       if (is.numeric(seed) && length(seed) == 1 && seed == as.integer(seed)) {
                         seed <<- seed
                       } else {
                         stop("Invalid seed")
                       }
                     },



                     # Generate sequence of hidden states
                     generate_sequence = function(length) {
                       set.seed(seed)
                       "Generate a sequence of hidden states of specified length"
                       if (length < 1) {
                         stop("Sequence length must be positive")
                       }

                       # Initialize sequence
                       sequence <- numeric(length)

                       # Generate first state based on initial probabilities
                       sequence[1] <- sample(1:n_states, 1, prob = init_prob)

                       # Generate subsequent states based on transition matrix
                       if (length > 1) {
                         for (t in 2:length) {
                           current_state <- sequence[t-1]
                           sequence[t] <- sample(1:n_states, 1, prob = trans_mat[current_state,])
                         }
                       }


                       # Convert numeric states to state names
                       named_sequence <- state_names[sequence]
                       return(list(
                         numeric_sequence = sequence,
                         named_sequence = named_sequence
                       ))
                     },

                     # Print HMM parameters
                     show = function() {
                       cat("Hidden Markov Model with", n_states, "states\n")
                       cat("\nState Names:\n")
                       print(state_names)
                       cat("\nTransition Matrix:\n")
                       print(trans_mat)
                       cat("\nInitial Probabilities:\n")
                       print(init_prob)
                     }
                   )
)

#' Title
#' HMM Gaussian Model Class
#'
#' @param hmm_model   # the base Hidden Markov Model used
#' @param state_variances # a vector of length K recording the variances of the Gaussian emission distribution
#' of each latent state
#' @param n_states  # a integer of number of latent states
#'
#' @return
#' This is a class of HMM Gaussian models to store the features of a Hidden Markov Model with Gaussian distribution as emission distribution
#' for continuous variables.
#' with given parameters.
#' @export
#'
#' @examples
#'

HMM_Gaussian_Model <- setRefClass("HMM_Gaussian_Model",
                                  fields = list(
                                    hmm_model = "HMM",            # Base HMM model
                                    state_variances = "numeric",   # Variances for each state's Gaussian emissions
                                    n_states = "numeric",          # Number of states (copied from HMM model)
                                    seed = "numeric"
                                  ),

                                  methods = list(

                                    # Initialize the Gaussian HMM
                                    initialize = function(hmm_model, state_variances) {
                                      "Initialize Gaussian HMM with base HMM model and state variances"

                                      # Validate input
                                      if (!is(hmm_model, "HMM")) {
                                        stop("hmm_model must be an instance of HMM class")
                                      }

                                      if (length(state_variances) != hmm_model$n_states) {
                                        stop("Number of variances must match number of states in HMM model")
                                      }

                                      if (any(state_variances <= 0)) {
                                        stop("All variances must be positive")
                                      }

                                      # Store parameters
                                      hmm_model <<- hmm_model
                                      state_variances <<- state_variances
                                      n_states <<- hmm_model$n_states
                                    },


                                    # Generate Gaussian observations
                                    generate_gaussian_observations = function(length) {
                                      "Generate sequence of Gaussian observations based on HMM states"

                                      # Generate latent sequence using base HMM
                                      latent_seq <- hmm_model$generate_sequence(length)

                                      # Generate observations
                                      observations <- numeric(length)

                                      for (i in 1:length) {
                                        current_state <- latent_seq$numeric_sequence[i]
                                        observations[i] <- rnorm(1,
                                                                 mean = current_state,
                                                                 sd = sqrt(state_variances[current_state]))
                                      }

                                      return(list(
                                        observations = observations,
                                        true_states = latent_seq$numeric_sequence,
                                        state_names = latent_seq$named_sequence
                                      ))
                                    },


                                    # Calculate emission probability density for a single observation
                                    emission_prob = function(observation, state) {
                                      "Calculate Gaussian emission probability density for given observation and state"
                                      dnorm(observation, mean = state, sd = sqrt(state_variances[state]))
                                    },


                                    # Calculate log-likelihood of a sequence
                                    calculate_log_likelihood = function(observations) {
                                      "Calculate log-likelihood of observation sequence using forward algorithm"

                                      T <- length(observations)
                                      alpha <- matrix(0, nrow = T, ncol = n_states)
                                      # Initialize forward probabilities (alpha) for t=1
                                      for (j in 1:n_states) {
                                        alpha[1, j] <- log(hmm_model$init_prob[j]) +
                                          log(emission_prob(observations[1], j))
                                      }
                                      # Forward algorithm recursion
                                      for (t in 2:T) {
                                        for (j in 1:n_states) {
                                          # Calculate alpha for each state
                                          logsum <- numeric(n_states)
                                          for (i in 1:n_states) {
                                            logsum[i] <- alpha[t-1, i] +
                                              log(hmm_model$trans_mat[i, j])
                                          }
                                          alpha[t, j] <- log(emission_prob(observations[t], j)) +
                                            max(logsum) +
                                            log(sum(exp(logsum - max(logsum))))
                                        }
                                      }


                                      # Calculate final log-likelihood
                                      logsum <- alpha[T, ]
                                      log_likelihood <- max(logsum) + log(sum(exp(logsum - max(logsum))))

                                      return(list(
                                        log_likelihood = log_likelihood,
                                        forward_probs = alpha
                                      ))
                                    },


                                    # Show method for displaying model parameters
                                    show = function() {
                                      cat("Gaussian Hidden Markov Model\n")
                                      cat("\nBase HMM Parameters:\n")
                                      hmm_model$show()
                                      cat("\nGaussian State Variances:\n")
                                      print(state_variances)
                                    }
                                  )
)

hmm_model = HMM(10)
