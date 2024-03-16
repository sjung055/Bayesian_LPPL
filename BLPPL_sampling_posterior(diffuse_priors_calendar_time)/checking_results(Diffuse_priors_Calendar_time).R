# <Paper>
# A Bayesian analysis of log-periodic precursors to financial crashes
# after sampling from posterior distribution

# A_{lp, c}: diffuse priors & calendar time model


title = "diffuse priors & calendar time"

## posterior mean values from the paper (A_{lp, c})
mu.paper = 0.000294
tau.paper = 15757
B.paper = 0.007195
C.paper = 0.498202
beta.paper = 0.361958
omega.paper = 6.429387
phi.paper = 3.124951
tc.paper = as.Date("1988-02-09")


# checking results
posterior_samples = list()
posterior_samples$mu = sample.mu
posterior_samples$tau = sample.tau
posterior_samples$B = sample.B
posterior_samples$C = sample.C
posterior_samples$beta = sample.beta
posterior_samples$omega = sample.omega
posterior_samples$phi = sample.phi
posterior_samples$tctN = sample.tctN
posterior_samples.names = names(posterior_samples)


## (1) acceptance rate
for(i in 1:length(posterior_samples)){
  acc_rate = length(unique(posterior_samples[[i]]))/iter_num
  print(sprintf("Acceptance Rate of %s: %f",
                posterior_samples.names[i], acc_rate))
}

## (2) ts plot
for(i in 1:length(posterior_samples)){
  ts.plot(posterior_samples[[i]],
          ylab = sprintf("%s", posterior_samples.names[i]), 
          main = sprintf("trace plot of %s (%s) before",
                         posterior_samples.names[i], title))
}

## (3) ACF plot
for(i in 1:length(posterior_samples)){
  acf(posterior_samples[[i]],
      main = sprintf("ACF plot of %s (%s) before",
                     posterior_samples.names[i], title))
}

## (4) burn-in & thinning
### burn-in
burn_in_period = 10000
posterior_samples$mu = sample.mu[(burn_in_period+1):iter_num]
posterior_samples$tau = sample.tau[(burn_in_period+1):iter_num]
posterior_samples$B = sample.B[(burn_in_period+1):iter_num]
posterior_samples$C = sample.C[(burn_in_period+1):iter_num]
posterior_samples$beta = sample.beta[(burn_in_period+1):iter_num]
posterior_samples$omega = sample.omega[(burn_in_period+1):iter_num]
posterior_samples$phi = sample.phi[(burn_in_period+1):iter_num]
posterior_samples$tctN = sample.tctN[(burn_in_period+1):iter_num]
### thinning
idx = seq(from = 1, to = (iter_num - burn_in_period), by = 3)
posterior_samples$mu = posterior_samples$mu[idx]
posterior_samples$tau = posterior_samples$tau[idx]
posterior_samples$B = posterior_samples$B[idx]
posterior_samples$C = posterior_samples$C[idx]
posterior_samples$beta = posterior_samples$beta[idx]
posterior_samples$omega = posterior_samples$omega[idx]
posterior_samples$phi = posterior_samples$phi[idx]
posterior_samples$tctN = posterior_samples$tctN[idx]


## (5) ts plot
for(i in 1:length(posterior_samples)){
  ts.plot(posterior_samples[[i]],
          ylab = sprintf("%s", posterior_samples.names[i]), 
          main = sprintf("trace plot of %s (%s) after",
                         posterior_samples.names[i], title))
}

## (6) ACF plot
for(i in 1:length(posterior_samples)){
  acf(posterior_samples[[i]],
      main = sprintf("ACF plot of %s (%s) after",
                     posterior_samples.names[i], title))
}

## (7) posterior mean & standard error
print("Posterior Mean & Standard Error")
for(i in 1:length(posterior_samples)){
  print(sprintf("<%s> Posterior Mean: %f", 
                posterior_samples.names[i], mean(posterior_samples[[i]])))
  print( sd(posterior_samples[[i]])/length(posterior_samples[[i]]) )
}
### t_c
tN + round(mean(posterior_samples$tctN))

## (8) histogram
hist(posterior_samples$mu, breaks = 50, freq = F,
     main = sprintf("mu (%s)", title), xlab = "mu",
     xlim = c(-0.004, 0.005))
abline(v = mu.paper, lwd = 2, col = "red")

hist(posterior_samples$tau, breaks = 50, freq = F,
     main = sprintf("tau (%s)", title), xlab = "tau",
     xlim = c(10000,20000))
abline(v = tau.paper, lwd = 2, col = "red")

hist(log(posterior_samples$B), breaks = 50, freq = F,
     main = sprintf("log(B) (%s)", title), xlab = "log(B)",
     xlim = c(-25,10))
abline(v = log(B.paper), lwd = 2, col = "red")

hist(posterior_samples$C, breaks = 50, freq = F,
     main = sprintf("C (%s)", title), xlab = "C",
     xlim = c(0,1))
abline(v = C.paper, lwd = 2, col = "red")

hist(posterior_samples$beta, breaks = 50, freq = F,
     main = sprintf("beta (%s)", title), xlab = "beta",
     xlim = c(0,1))
abline(v = beta.paper, lwd = 2, col = "red")

hist(posterior_samples$omega, breaks = 50, freq = F,
     main = sprintf("omega (%s)", title), xlab = "omega",
     xlim = c(0,20))
abline(v = omega.paper, lwd = 2, col = "red")
     
hist(posterior_samples$phi, breaks = 50, freq = F,
     main = sprintf("phi (%s)", title), xlab = "phi",
     xlim = c(0,2*pi))
abline(v = phi.paper, lwd = 2, col = "red")
     
hist(posterior_samples$tctN, breaks = 50, freq = F,
     main = sprintf("tc - tN (%s)", title), xlab = "tc - tN",
     xlim = c(0,1000))
abline(v = tctN.paper, lwd = 2, col = "red")

