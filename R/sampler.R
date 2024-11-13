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

stan_fit <- function(model_path, stan_data){
  hmm_model <- stan_model(model_code = readLines(model_path))
  fit <- sampling(hmm_model,
                  data = stan_data,
                  chains = 4,
                  iter = 2000,
                  warmup = 1000,
                  init = rep(list(init_fun()), 4),  # Use our initialization function
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  return(fit)
}

stan_to_model <- function(model_path){

}
