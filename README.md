# Bayesian Log Periodic Power Law(LPPL) Model

## Description
This repository contains the implementation of algorithms outlined in the paper ([George Chang & James Feigenbaum (2006) A Bayesian analysis of log-periodic precursors to financial crashes, Quantitative Finance, 6:1, 15-36](https://doi.org/10.1080/14697680500511017)). Since the paper was published in 2006, the algorithms it proposes have not been available in a coded form. In response, I developed code to replicate the paper’s results. My methodology relied on the Markov Chain Monte Carlo (MCMC) method, specifically utilizing the Metropolis-Hastings Algorithm within Gibbs Sampler (referred to as <ins>**Metropolis within Gibbs**</ins>), to accurately reproduce the findings of the original study.

- data: the S&P 500 index data spanning from 01/02/1980 to 12/30/1988 ([`S&P_500_Index_01021980_12301988_preprocessed.csv`](data/S&P_500_Index_01021980_12301988_preprocessed.csv))

(For comprehensive details and structures about Metropolis within Gibbs employed in this code, please refer to [Implementation of “A Bayesian analysis of log-periodic precursors to financial crashes” paper in R](./Implementation%20of%20“A%20Bayesian%20analysis%20of%20log-periodic%20precursors%20to%20financial%20crashes”%20paper%20in%20R.pdf))


## Usage
1. Execute [`utils_BLPPL_posterior.R`](./utils_BLPPL_posterior.R).
1. Select a combination of prior distributions and time models.
    - prior distributions: diffuse priors, tight priors
    - time models: calendar time model, market time model
1. Navigate to **the folder corresponding to your chosen combination**.  
    ex. For diffuse priors & calendar time model, go to BLPPL_sampling_posterior(diffuse_priors_calendar_time).
1. Execute `BLPPL_posterior(chosen combination).R` to obtain posterior samples for each parameter.  
    ex. For diffuse priors & calendar time model, run BLPPL_posterior(Diffuse_priors_Calendar_time).R.
1. To analyze the posterior distribution, execute `checking_results(chosen combination).R`.  
    ex. For diffuse priors & calendar time model, run checking_results(Diffuse_priors_Calendar_time).R.
1. A posterior sample of parameters at each iteration allows for plotting a line using the LPPL Model equation. To address uncertainty, the Bayesian approach advocates using a credible interval. I have thus developed a code ([`posterior_credible_interval.R`](./posterior_credible_interval.R)) that plots a 95% credible interval alongside the actual S&P 500 index data. This feature represents a novel contribution of my code, as the original paper did not discuss the use of a 95% credible interval.

## Contact
For any inquiries, you can reach me at [sjung055@yonsei.ac.kr](sjung055@yonsei.ac.kr).
