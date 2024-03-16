# <Paper>
# A Bayesian analysis of log-periodic precursors to financial crashes
# posterior distribution part

# B_{lp, m}: tight priors & market time model (not in the paper)

# 1. initial values & prior setting
# Black Monday: 10/19/1987
## initial values(they are set to NLLS estimates except for mu and tau)
mu.init = 0
tau.init = 15000
B.init = 0.0130
C.init = 0.966
beta.init = 0.580
omega.init = 5.711
phi.init = 4.845
# tc.init = as.Date("1987-10-20")

## posterior mean values from the paper (B_{lp, c})
mu.paper = 0.000032
tau.paper = 15698
B.paper = 0.012553
C.paper = 0.721244
beta.paper = 0.530319
omega.paper = 5.940680
phi.paper = 3.818973
# tc.paper = as.Date("1987-10-20")

## parameters of tight priors
mu.mu = 0.0
mu.sig2 = 10^(-6)
tau.shape = 1.0
tau.rate = 10^(-5)
B.shape = 100.0
B.rate = 7613.8
C.min = 0
C.max = 1
beta.shape1 = 41.3834
beta.shape2 = 29.9228
omega.shape = 16.0
omega.rate = 2.5
phi.min = 0
phi.max = 2*pi
tctN.shape = 100.0
tctN.rate = 25.0

########################################################################

# 2. data setting
## S&P 500 index data for tight priors: 1/3/1983 ~ 10/16/1987
tight_priors.data = snp_data[("1983-01-03" <= snp_data$Date) &
                               (snp_data$Date <= "1987-10-16"), ]
# N+1 number of days
#   (t_0: dtight_priors.data[1], t_2: tight_priors.data[2], ...,
#    t_N: tight_priors.data[N+1])
N1 = dim(tight_priors.data)[1]

## taking log to the price to make q(t)
### q(t_{i+1}): [q(t_1), q(t_2), ..., q(t_{N})]
q_ti1 = log(tight_priors.data[2:N1, ]$"Close")
### q(t_i): [q(t_0), q(t_1), ..., q(t_{N-1})]
q_ti = log(tight_priors.data[1:(N1-1), ]$"Close")

## market time model
### t_{i+1}: [t_1, t_2, ..., t_N]
ti1 = seq(from = 1, to = N1-1, by = 1)
### t_i: [t_0, t_2, ..., t_{N-1}]
ti = seq(from = 0, to = N1-2, by = 1)

### t
tN = calc_market_time(snp_data, as.Date("1983-01-03"), as.Date("1987-10-16"))
tc.init = calc_market_time(snp_data, as.Date("1983-01-03"),
                           as.Date("1987-10-20"))
tc.paper = calc_market_time(snp_data, as.Date("1983-01-03"),
                           as.Date("1987-10-20"))
tctN.init = tc.init - tN # t_c - t_N
tctN.paper = tc.paper - tN

########################################################################

# 3. hyper parameter: step size of proposal distribution

# sigma value for proposal distribution N(mu^t, p.mu.step_size^2)
p.mu.step_size = 0.00075
# sigma value for proposal distribution N(tau^t, p.tau.step_size^2)
p.tau.step_size = 2500
# sigma value for proposal distribution N(B^t, p.B.step_size^2)
p.B.step_size = 0.005
# sigma value for proposal distribution N(C^t, p.C.step_size^2)
p.C.step_size = 0.75
# sigma value for proposal distribution N(beta^t, p.beta.step_size^2)
p.beta.step_size = 0.15
# sigma value for proposal distribution N(omega^t, p.omega.step_size^2)
p.omega.step_size = 1
# sigma value for proposal distribution N(phi^t, p.phi.step_size^2)
p.phi.step_size = 2.5
# sigma value for proposal distribution N(tctN^t, p.tctN.step_size^2)
p.tctN.step_size = 1.5


iter_num = 100000 # total number of iterations

sample.mu = rep(NA, iter_num)
sample.tau = rep(NA, iter_num)
sample.B = rep(NA, iter_num)
sample.C = rep(NA, iter_num)
sample.beta = rep(NA, iter_num)
sample.omega = rep(NA, iter_num)
sample.phi = rep(NA, iter_num)
sample.tctN = rep(NA, iter_num)
sample.mu[1] = mu.init
sample.tau[1] = tau.init
sample.B[1] = B.init
sample.C[1] = C.init
sample.beta[1] = beta.init
sample.omega[1] = omega.init
sample.phi[1] = phi.init
sample.tctN[1] = tctN.init

########################################################################
# 4. implementing "MH algorithm within Gibbs sampler"
start_time <- Sys.time()
for(t in 1:(iter_num-1)){
  mu.t = sample.mu[t]
  tau.t = sample.tau[t]
  B.t = sample.B[t]
  C.t = sample.C[t]
  beta.t = sample.beta[t]
  omega.t = sample.omega[t]
  phi.t = sample.phi[t]
  tctN.t = sample.tctN[t]
  
  tc.t = tctN.t + tN
  
  
  # sampling mu
  p.mu = rnorm(n = 1, mean = mu.t, sd = p.mu.step_size) # proposed mu
  p.param = c(p.mu, tau.t, B.t, C.t, beta.t, omega.t, phi.t)
  param.t = c(mu.t, tau.t, B.t, C.t, beta.t, omega.t, phi.t)
  log_numerator = dnorm(x = p.mu,
                        mean = mu.mu, sd = sqrt(mu.sig2), log = T) +
                  sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                              p.param, tc.t))
  log_denominator = dnorm(x = mu.t,
                          mean = mu.mu, sd = sqrt(mu.sig2), log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                param.t, tc.t))
  log_alpha = log_numerator - log_denominator
  mu.t1 = ifelse(log_alpha>log(runif(n=1)), p.mu, mu.t)
  sample.mu[t+1] = mu.t1
  
  # sampling tau
  p.tau = rnorm(n = 1, mean = tau.t, sd = p.tau.step_size) # proposed tau
  p.param = c(mu.t1, p.tau, B.t, C.t, beta.t, omega.t, phi.t)
  param.t = c(mu.t1, tau.t, B.t, C.t, beta.t, omega.t, phi.t)
  if(p.tau<=0){
    # reject
    sample.tau[t+1] = tau.t
    tau.t1 = tau.t # tau_{t+i}
  }else{
    # p.tau > 0
    log_numerator = dgamma(x = p.tau,
                           shape = tau.shape, rate = tau.rate, log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                p.param, tc.t))
    log_denominator = dgamma(x = tau.t,
                             shape = tau.shape, rate = tau.rate, log = T) +
                      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    tau.t1 = ifelse(log_alpha>log(runif(n=1)), p.tau, tau.t)
    sample.tau[t+1] = tau.t1
  }
  
  # sampling B
  p.B = rnorm(n = 1, mean = B.t, sd = p.B.step_size) # proposed B
  p.param = c(mu.t1, tau.t1, p.B, C.t, beta.t, omega.t, phi.t)
  param.t = c(mu.t1, tau.t1, B.t, C.t, beta.t, omega.t, phi.t)
  if(p.B<=0){
    # reject
    sample.B[t+1] = B.t
    B.t1 = B.t # B_{t+i}
  }else{
    # p.B > 0
    log_numerator = dgamma(x = p.B,
                           shape = B.shape, rate = B.rate, log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                p.param, tc.t))
    log_denominator = dgamma(x = B.t,
                             shape = B.shape, rate = B.rate, log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    B.t1 = ifelse(log_alpha>log(runif(n=1)), p.B, B.t)
    sample.B[t+1] = B.t1
  }
  
  # sampling C
  p.C = rnorm(n = 1, mean = C.t, sd = p.C.step_size) # proposed C
  p.param = c(mu.t1, tau.t1, B.t1, p.C, beta.t, omega.t, phi.t)
  param.t = c(mu.t1, tau.t1, B.t1, C.t, beta.t, omega.t, phi.t)
  if( (p.C<0) | (1<p.C) ){
    # reject
    sample.C[t+1] = C.t
    C.t1 = C.t # C_{t+i}
  }else{
    # 0 <= p.C <= 1
    log_numerator = dunif(x = p.C,
                          min = C.min, max = C.max, log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                p.param, tc.t))
    log_denominator = dunif(x = C.t,
                            min = C.min, max = C.max, log = T) +
                      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    C.t1 = ifelse(log_alpha>log(runif(n=1)), p.C, C.t)
    sample.C[t+1] = C.t1
  }
  
  # sampling beta
  p.beta = rnorm(n = 1, mean = beta.t, sd = p.beta.step_size) # proposed beta
  p.param = c(mu.t1, tau.t1, B.t1, C.t1, p.beta, omega.t, phi.t)
  param.t = c(mu.t1, tau.t1, B.t1, C.t1, beta.t, omega.t, phi.t)
  if( (p.beta<0) | (1<p.beta) ){
    # reject
    sample.beta[t+1] = beta.t
    beta.t1 = beta.t # beta_{t+i}
  }else{
    # 0 <= p.beta <= 1
    log_numerator = dbeta(x = p.beta, shape1 = beta.shape1,
                          shape2 = beta.shape2, log = T) +
                    sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                p.param, tc.t))
    log_denominator = dbeta(x = beta.t, shape1 = beta.shape1,
                            shape2 = beta.shape2, log = T) +
                      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    beta.t1 = ifelse(log_alpha>log(runif(n=1)), p.beta, beta.t)
    sample.beta[t+1] = beta.t1
  }
  
  # sampling omega
  p.omega = rnorm(n = 1, mean = omega.t, sd = p.omega.step_size) # proposed omega
  p.param = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, p.omega, phi.t)
  param.t = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, omega.t, phi.t)
  if(p.omega<=0){
    # reject
    sample.omega[t+1] = omega.t
    omega.t1 = omega.t # omega_{t+i}
  }else{
    # p.omega > 0
    log_numerator = dgamma(x = p.omega, shape = omega.shape,
                           rate = omega.rate, log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  p.param, tc.t))
    log_denominator = dgamma(x = omega.t, shape = omega.shape,
                             rate = omega.rate, log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    omega.t1 = ifelse(log_alpha>log(runif(n=1)), p.omega, omega.t)
    sample.omega[t+1] = omega.t1
  }
  
  # sampling phi
  p.phi = rnorm(n = 1, mean = phi.t, sd = p.phi.step_size) # proposed phi
  p.param = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, omega.t1, p.phi)
  param.t = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, omega.t1, phi.t)
  if( (p.phi<phi.min) | (phi.max<p.phi) ){
    # reject
    sample.phi[t+1] = phi.t
    phi.t1 = phi.t # phi_{t+i}
  }else{
    # 0 <= p.phi <= 2*pi
    log_numerator = dunif(x = p.phi,
                          min = phi.min, max = phi.max, log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  p.param, tc.t))
    log_denominator = dunif(x = phi.t,
                            min = phi.min, max = phi.max, log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    phi.t1 = ifelse(log_alpha>log(runif(n=1)), p.phi, phi.t)
    sample.phi[t+1] = phi.t1
  }
  
  # sampling tctN
  p.tctN = rnorm(n = 1, mean = tctN.t, sd = p.tctN.step_size) # proposed tctN
  p.tc = p.tctN + tN
  
  p.param = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, omega.t1, phi.t1)
  param.t = c(mu.t1, tau.t1, B.t1, C.t1, beta.t1, omega.t1, phi.t1)
  if( p.tc < tN ){
    # reject
    sample.tctN[t+1] = tctN.t
    tctN.t1 = tctN.t
  }else{
    # tN <= p.tc
    log_numerator = dgamma(x = p.tctN, shape = tctN.shape, rate = tctN.rate,
                           log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  p.param, p.tc))
    log_denominator = dgamma(x = tctN.t, shape = tctN.shape, rate = tctN.rate,
                             log = T) +
      sum(log_probability_density(q_ti1, q_ti, ti1, ti,
                                  param.t, tc.t))
    log_alpha = log_numerator - log_denominator
    tctN.t1 = ifelse(log_alpha>log(runif(n=1)), p.tctN, tctN.t)
    sample.tctN[t+1] = tctN.t1
  }
}
end_time <- Sys.time()
end_time - start_time
# end of MCMC
########################################################################

