# Drawing plot with posterior samples

# At first, load posterior sample data

model_name = "Tight priors & Calendar time"

calendar_days= snp_data$Date[snp_data$Date <= tN ]

#########################################################################
# <Draw a line using posterior mean>

A = 5.7

log_snp_500_index_col = "black"
post_mean_line_col = "red"
log_snp_500_index_type = "l"
log_snp_500_index_lty = "solid"
post_mean_line_type = "l"
post_mean_line_lty = "dashed"

plot(calendar_days, log(snp_data$Close[snp_data$Date <= tN ]),
     type = log_snp_500_index_type, lty = log_snp_500_index_lty,
     col = log_snp_500_index_col,
     xlab = "Time", ylab = "q(t)", main = model_name)
lines(calendar_days,
      A + H_t(mean(posterior_samples$B), mean(posterior_samples$C),
          mean(posterior_samples$beta), mean(posterior_samples$omega),
          mean(posterior_samples$phi),
          as.numeric(tN + round(mean(posterior_samples$tctN)) - calendar_days)),
      col = post_mean_line_col,
      type = post_mean_line_type, lty = post_mean_line_lty)
legend("topleft",
       legend=c("log(S&P 500 index)", "line with posterior mean valeus"),
       lty = c(log_snp_500_index_lty, post_mean_line_lty),
       col = c(log_snp_500_index_col, post_mean_line_col),
       border = "white", box.lty = 1, cex = 1)


#########################################################################

# <95% credible interval>

A = 5.7

# H_t.posterior is a 30,000*1970 matrix, and its i-th row stores
# from H_1 to H_{t_N} generated with i-th posterior sample parameters.
H_t.posterior = matrix(NA, nrow = length(posterior_samples[[1]]),
                      ncol = length(calendar_days))
for( t in seq(length(posterior_samples[[1]])) ){
  H_t.posterior[t, ] = H_t(posterior_samples$B[t], posterior_samples$C[t],
                          posterior_samples$beta[t], posterior_samples$omega[t],
                          posterior_samples$phi[t],
                          as.numeric(tN + posterior_samples$tctN[t] - calendar_days))
}

## getting 2.5 percentile and 97.5 percentile value for each day
credible_interval = apply(H_t.posterior, 2,
                          function(column) quantile(column, probs = c(0.025, 0.975)))

## drawing log(price) graph with 2.5 and 9.75 credible interval graph
log_snp_500_index_col = "black"
credible_interval_col = "red"
log_snp_500_index_type = "l"
log_snp_500_index_lty = "solid"
credible_interval_type = "l"
credible_interval_lty = "dashed"

plot(calendar_days, log(snp_data$Close[snp_data$Date <= tN ]),
     type = log_snp_500_index_type, lty = log_snp_500_index_lty,
     ylim = c(4,5.8), xlab = "Time", ylab = "q(t)",
     main = model_name, col = log_snp_500_index_col)
# lower interval
lines(calendar_days, A + credible_interval[1, ],
      col = credible_interval_col,
      type = 'l', lty = credible_interval_lty)
# upper interval
lines(calendar_days, A + credible_interval[2, ],
      col = credible_interval_col,
      type = credible_interval_type,
      lty = credible_interval_lty)
legend("bottomright",
       legend=c("log(S&P 500 index)", "95% credible interval"),
       lty = c(log_snp_500_index_lty, credible_interval_lty),
       col = c(log_snp_500_index_col, credible_interval_col),
       border = "white", box.lty = 1, cex = 0.8)




