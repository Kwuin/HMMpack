# Hidden Markov Model for sequential data
This is a R package for implementing Hidden Markov Models(HMM) for sequential data. It includes functionalities of generating synthetic data with arbitray transition matrices, initial distribution and differnet emission distribution for continuous data. It can also be used to do Bayesian inference using Rstan to do Hamiltonian Monte Carlo sampling method. Calculation of log likelihood for given sequential data and model is provided. 

## Installation
Run the following code in Rstudio: 
```{r}
devtools::install_github("Kwuin/HMMpack")
```

## Remaning part
1. Result visualization: the sampling result from Rstan should be organized and visualized properly. It includes the MAP point estimation of initial distribution, transition matrix and also an uncertainty quantification report of each parameter. 
2. Cross validation: the main goal of modeling sequential data is to do the predictions. We would implement a cross validation with a sliding-window based iterative inference and a one shot prediction. Then we would show the error between the prediction and the real data. The error can be $L_0$, $L_1$, $L_2$, cross-entropy, etc. A plot report would be generated. 
3. On condition on sufficiency of time, we may explore more inference methods like expectation maximization and variational inference. 