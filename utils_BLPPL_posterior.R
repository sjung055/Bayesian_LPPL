# <Paper>
# A Bayesian analysis of log-periodic precursors to financial crashes

# utility functions for importance sampling part

rm(list = ls())
set.seed(123)

## data
library(tidyverse)
snp_data = read_csv("./data/S&P_500_Index_01021980_12301988_preprocessed.csv")
snp_data = snp_data[order(snp_data$Date, decreasing=FALSE),]

### check the data
head(snp_data)

plot(x = snp_data$Date, y = snp_data$Close, type = 'l',
     xlab = "Date", ylab = "S&P 500 index",
     main = sprintf("S&P 500 (%s ~ %s)",
                    snp_data$Date[1],
                    snp_data$Date[dim(snp_data)[1]]))


# Function for calculating market time between two days
## ex.
## 10.19.1987(Mon) - 10.16.1987(Fri) = 1
## 10.16.1987(Fri) - 10.14.1987(Wed) = 2
calc_market_time = function(snp_data, start_date, end_date){
  # start_date, end_date: Date type
  bewteen_days = snp_data[(start_date <= snp_data$Date) &
                          (snp_data$Date <= end_date), ]
  return(dim(bewteen_days)[1] - 1)
}



# Functions from the paper

## (1) function for calculating H(t)
H_t = function(B, C, beta, omega, phi, tct){
  # tct: t_c - t
  part1 = -B * tct^beta
  part2 = 1 + C/sqrt(1+(omega/beta)^2) * cos(omega * log(tct) + phi)
  return(part1*part2) # H(t)
}


## (2) log[ p( q_{t_{i+1}}|q_{t_i}, \theta ) ]
log_probability_density = function(q_ti1, q_ti, ti1, ti, params, tc){
  mu =  params[1]; tau = params[2]; B = params[3]; C = params[4];
  beta = params[5]; omega = params[6]; phi = params[7]
  # q_ti1: q_{t_{i+1}}
  # q_ti: q_{t_i}
  # ti1: t_{i+1}
  # ti: t_i
  # tct_i1: t_c - t_{i+1}
  # tct_i: t_c - t_i
  tct_i1 = tc - ti1
  tct_i = tc - ti
  
  part1 = (1/2)*log(tau)
  part2 = -(1/2)*log(2*pi*as.numeric(ti1-ti))

  H_ti1 = H_t(B, C, beta, omega, phi, as.numeric(tct_i1)) # H(t_{i+1})
  H_ti = H_t(B, C, beta, omega, phi, as.numeric(tct_i)) # H(t_i)
  # dH: Delta H(t_i, t_{i+1}; xi)
  dH = H_ti1 - H_ti

  part3 = -tau*(q_ti1-q_ti -mu*as.numeric(ti1-ti) - dH)^2 /
    (2*as.numeric(ti1-ti))

  return(part1+part2+part3)
}


## (3) log crash probabilities
log_crash_prob = function(t0_to_t_N1, kappa, params, tc){
  mu =  params[1]; tau = params[2]; B = params[3]; C = params[4];
  beta = params[5]; omega = params[6]; phi = params[7]
  # t0_to_t_N1: date from t_0 to t_{N+1}
  # ex. t_0 = t0_to_t_N1[1], t_1 = t0_to_t_N1[2], ..., t_N = t0_to_t_N1[N+1],
  #     t_{N+1} = t0_to_t_N1[N+2]
  N = length(t0_to_t_N1) - 2
  
  # t1_to_t_N: date from t_1 to t_N
  # t0_to_tN_1: date from t_0 to t_{N-1}
  t1_to_t_N = t0_to_t_N1[2:(N+1)]
  t0_to_tN_1 = t0_to_t_N1[1:N]
  # tct_i: t_c - t_i
  # tct_i_1: t_c - t_{i-1}
  tct_i = tc - t1_to_t_N
  tct_i_1 = tc - t0_to_tN_1

  H_ti = H_t(B, C, beta, omega, phi, as.numeric(tct_i)) # H(t_i)
  H_ti_1 = H_t(B, C, beta, omega, phi, as.numeric(tct_i_1)) # H(t_{i-1})
  # dH: Delta H(t_{i-1}, t_i; xi)
  dH = H_ti - H_ti_1
  
  part1 = -sum((1/kappa)*dH)
  
  tct_min = tc - min(t0_to_t_N1[N+2], tc)
  H_min_tN1_tc = H_t(B, C, beta, omega, phi, as.numeric(tct_min)) # H(min{t_{N+1}, t_c})
  tctN = tc - t0_to_t_N1[N+1] # t_c - t_N
  H_tN = H_t(B, C, beta, omega, phi, as.numeric(tctN)) # H(t_N)
  # dH: Delta H(t_N, min{t_{N+1}, t_c}; xi)
  dH = H_min_tN1_tc - H_tN
  part2 = log(1 - exp(-(1/kappa) * dH))
  ########################################################
  ifelse(1 - exp(-(1/kappa) * dH) <= 0, return(-Inf), return(part1+part2))
  ########################################################
  # return(part1+part2)
}

