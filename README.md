# Bayesian Log Periodic Power Law(LPPL) Model

## Description
This code is for reproducing the sampled posterior distribution results from the paper ([George Chang & James Feigenbaum (2006) A Bayesian analysis of log-periodic precursors to financial crashes, Quantitative Finance, 6:1, 15-36](https://doi.org/10.1080/14697680500511017)).

(For detailed information about Metropolis-Hastings Algorithm within
the Gibbs Sampler which is used in this code, please refer to [Implementation of “A Bayesian analysis of log-periodic precursors to financial crashes” paper in R](./Implementation%20of%20“A%20Bayesian%20analysis%20of%20log-periodic%20precursors%20to%20financial%20crashes”%20paper%20in%20R.pdf))


## Usage
1. Run [`utils_BLPPL_posterior.R`](./utils_BLPPL_posterior.R).
1. Choose a combination of prior distributions and time models.
    - prior distributions: diffuse priors, tight priors
    - time models: calendar time model, market time model
1. Go to **the folder of the combination**.  
    ex. For diffuse priors & calendar time model, go to BLPPL_sampling_posterior(diffuse_priors_calendar_time).
1. Run `BLPPL_posterior(chosen combination).R`. Then, you will get posterior samples of each parameter.  
    ex. For diffuse priors & calendar time model, run BLPPL_posterior(Diffuse_priors_Calendar_time).R.
1. To check the posterior distribution, run `checking_results(chosen combination).R`.  
    ex. For diffuse priors & calendar time model, run checking_results(Diffuse_priors_Calendar_time).R.
