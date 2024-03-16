# <Paper>
# A Bayesian analysis of log-periodic precursors to financial crashes
# after sampling from posterior distribution

# A_{lp, m}: diffuse priors & market time model

title = "diffuse priors & market time"

## posterior mean values from the paper (A_{lp, m})
mu.paper = 0.000420
tau.paper = 13437
B.paper = 0.007359
C.paper = 0.500030
beta.paper = 0.384268
omega.paper = 6.425584
phi.paper = 3.149531
# tc.paper = as.Date("1988-03-31")
### t
tN = calc_market_time(snp_data, as.Date("1983-01-03"), as.Date("1987-10-16"))
tc.init = calc_market_time(snp_data, as.Date("1983-01-03"),
                           as.Date("1987-10-20"))
tc.paper = calc_market_time(snp_data, as.Date("1983-01-03"),
                            as.Date("1988-03-31"))
tctN.init = tc.init - tN # t_c - t_N
tctN.paper = tc.paper - tN


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

## (7) posterior mean
print("Posterior Mean & Standard Error")
for(i in 1:length(posterior_samples)){
  print(sprintf("<%s> Posterior Mean: %f", 
                posterior_samples.names[i], mean(posterior_samples[[i]])))
  print( sd(posterior_samples[[i]])/length(posterior_samples[[i]]) )
}
### meaning of tc
tctN.post_mean = round(mean(posterior_samples$tctN))
end_date_idx = which("1987-10-16" == snp_data$Date)
snp_data$Date[end_date_idx+tctN.post_mean] # critical date


## (8) histogram
colors_vec = c("red", "dodgerblue")
lty_vec = c("solid", "dotted")
hist(posterior_samples$mu, breaks = 50, freq = F,
     main = sprintf("mu (%s)", title), xlab = "mu",
     xlim = c(-0.004, 0.005))
abline(v = mu.paper, lty = lty_vec[1], lwd = 2, col = colors_vec[1])
abline(v = mean(posterior_samples$mu), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = lty_vec, lwd = 2, cex = 0.7,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$tau, breaks = 50, freq = F,
     main = sprintf("tau (%s)", title), xlab = "tau",
     xlim = c(10000,20000))
abline(v = tau.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$tau), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.8,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(log(posterior_samples$B), breaks = 50, freq = F,
     main = sprintf("log(B) (%s)", title), xlab = "log(B)",
     xlim = c(-25,10))
abline(v = log(B.paper), lty = lty_vec[2], lwd = 2, col = "red")
abline(v = mean(log(posterior_samples$B)), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topleft", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.7,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$C, breaks = 50, freq = F,
     main = sprintf("C (%s)", title), xlab = "C",
     xlim = c(0,1))
abline(v = C.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$C), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.6,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$beta, breaks = 50, freq = F,
     main = sprintf("beta (%s)", title), xlab = "beta",
     xlim = c(0,1))
abline(v = beta.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$beta), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.8,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$omega, breaks = 50, freq = F,
     main = sprintf("omega (%s)", title), xlab = "omega",
     xlim = c(0,20))
abline(v = omega.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$omega), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.8,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$phi, breaks = 50, freq = F,
     main = sprintf("phi (%s)", title), xlab = "phi",
     xlim = c(0,2*pi))
abline(v = phi.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$phi), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("top", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.5,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)

hist(posterior_samples$tctN, breaks = 50, freq = F,
     main = sprintf("tc - tN (%s)", title), xlab = "tc - tN",
     xlim = c(0,1000))
abline(v = tctN.paper, lty = lty_vec[1], lwd = 2, col = "red")
abline(v = mean(posterior_samples$tctN), lty = lty_vec[2],
       lwd = 2, col = colors_vec[2])
legend("topright", col = colors_vec, border = "white", box.lty = 0,
       lty = "solid", lwd = 2, cex = 0.8,
       legend = c("posterior mean(paper)",
                  "posterior mean(this algorithm)"),
)
