#### Question 2
library(mvtnorm)
library(ivmodel)
library(ggplot2)

# We will run 200000 simulations
n_sims = 1000

# Specified parameters
beta = 0
delta = 0.5
n = 100
sig_mat = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)

# Set up output df (will parallelize later)
to_collect <- c("b_tsls", "b_liml", "var_x", "var_u", "rsq", "var_v", "cov_uv", "waa_qz", "waa_delta")
simulated_data <- data.frame(matrix(ncol = length(to_collect), nrow = n_sims))
colnames(simulated_data) <- to_collect

for (i in  1:n_sims) {
  if (((i - 1) %% 100 == 0) & (i != 1)) {
    cat("Finished iteration:", i - 1, "\n")
  } 
  
  # Generate data
  z_i = rnorm(n)
  uv = rmvnorm(n, sigma = sig_mat)
  x_i = (delta/(sqrt(n))) * z_i + uv[,2]
  y_i = beta * x_i + uv[,1]
  
  # Get tsls and liml estimates
  mod <- ivmodel(Y = y_i, D = x_i, Z = z_i)
  
  tsls_est <- coef(mod)["TSLS", "Estimate"]
  liml_est <- coef(mod)["LIML", "Estimate"]
  
  # For 2SLS normal asymptotics, need var(x), var(u), R^2(x~z)
  # var(x)
  var_x = var(x_i)
  
  # var(u)
  u_i_hat = y_i - tsls_est * x_i
  var_u_hat = var(u_i_hat)
  var_u = var(uv[,1])
  
  # R^2(x~z)
  rsq = summary(lm(x_i ~ z_i))$r.squared
  
  # For weak asymptotic approx, need var(u), var(v), cov(u, v), delta, qz
  # var(u) from above
  
  # QZ 
  Z_ = as.matrix(z_i)
  waa_qz = t(Z_) %*% Z_ / n
  
  # Delta is coef of regression of x~z multiplied by sqrt(n)
  
  waa_delta = summary(lm(x_i ~ z_i))$coef["z_i", "Estimate"] * sqrt(n)
  
  # var(v)
  v_i_hat = x_i - (waa_delta / sqrt(n)) * z_i
  var_v_hat = var(v_i_hat)
  var_v = var(uv[,2])
  
  # cov(u, v)
  cov_uv_hat = cov(u_i_hat, v_i_hat)
  cov_uv = cov(uv[,1], uv[,2])
  
  simulated_data[i,] = c(tsls_est, liml_est, var_x, var_u, rsq, var_v, cov_uv, waa_qz, waa_delta)
}

# Then the normal variance is this
simu_norm_var = mean(simulated_data$var_u) / (mean(simulated_data$var_x) * mean(simulated_data$rsq))
# But, we have to shrink by dividing by sqrt(n), which divides variance by n
std_asym_var = simu_norm_var / n

# Which we now simulate from the normal
simulated_data$normal_app = rnorm(n_sims, 0, sqrt(simu_norm_var/10))

# For the weak asymptotic approach, we estimate v-cov matrix and QZ*delta
sig_hat = matrix(c(mean(simulated_data$var_u), mean(simulated_data$cov_uv),
                   mean(simulated_data$cov_uv), mean(simulated_data$var_v)),
                 nrow = 2, ncol = 2)
qz_delt = mean(simulated_data$waa_qz) * mean(simulated_data$waa_delta)

# And we can simulate using the weak asymptotic distribution
uv_draws = rmvnorm(n_sims, sigma = sig_hat)

simulated_data$waa_app = uv_draws[,1] / (qz_delt + uv_draws[,2])

ggplot(simulated_data) + 
  geom_density(aes(x = b_tsls, color = "Exact Distribution")) + 
  geom_density(aes(x = normal_app, color = "Standard Asymptotics Normal Approximation")) + 
  geom_density(aes(x = waa_app, color = "Weak Instrument Asymptotics Approximation")) +
  scale_x_continuous(limits = c(-6, 6)) + theme(legend.position = "bottom")



