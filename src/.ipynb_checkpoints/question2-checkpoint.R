#### Question 2
library(mvtnorm)
library(ivmodel)
library(ggplot2)
library(doParallel)
library(foreach)
library(doRNG)

rm(list = ls())

# We will run 200000 simulations
n_sims = 200000

ncore <- detectCores()
cl <- makeCluster(ncore - 2, type = "PSOCK") #ncore - 1 usually
registerDoParallel(cl)

# From random.org
set.seed(859255299)
time_parallel <- system.time({
  simulation_data <- foreach(i = 1:n_sims
                             , .combine = 'rbind'
                             , .packages = c("ivmodel", "mvtnorm")
  ) %dorng% {
    
    # Specified parameters
    beta = 0
    delta = 0.5
    n = 100
    sig_mat = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    
    # Generate data
    z_i = rnorm(n)
    uv = rmvnorm(n, sigma = sig_mat)
    x_i = (delta/(sqrt(n))) * z_i + uv[,2]
    y_i = beta * x_i + uv[,1]
    
    # Get tsls and liml estimates
    mod <- ivmodel(Y = y_i, D = x_i, Z = z_i)
    
    b_tsls <- coef(mod)["TSLS", "Estimate"]
    b_liml <- coef(mod)["LIML", "Estimate"]
    
    # For 2SLS normal asymptotics, need var(x), var(u), R^2(x~z)
    # var(x)
    var_x = var(x_i)
    
    # var(u)
    u_i_hat = y_i - b_tsls * x_i
    var_u_hat = var(u_i_hat)
    var_u = var(uv[,1])
    
    # R^2(x~z)
    rsq = summary(lm(x_i ~ z_i))$r.squared
    
    # For weak asymptotic approx, need var(u), var(v), cov(u, v), delta, qz
    # var(u) from above
    
    # QZ 
    Z_ = as.matrix(z_i)
    wa_qz = t(Z_) %*% Z_ / n
    
    # Delta is coef of regression of x~z multiplied by sqrt(n)
    
    wa_delta = summary(lm(x_i ~ z_i))$coef["z_i", "Estimate"] * sqrt(n)
    
    # var(v)
    v_i_hat = x_i - (wa_delta / sqrt(n)) * z_i
    var_v_hat = var(v_i_hat)
    var_v = var(uv[,2])
    
    # cov(u, v)
    cov_uv_hat = cov(u_i_hat, v_i_hat)
    cov_uv = cov(uv[,1], uv[,2])
    
    # The standard error calculated when we run 2SLS should be based off the asymptotic model
    tsls_se <- coef(mod)["TSLS", "Std. Error"]
    
    simulated_data = c(b_tsls, b_liml, var_x, var_u, rsq, var_v, cov_uv, wa_qz, wa_delta, tsls_se)
    
    simulated_data
  }
})
stopCluster(cl)

cat("Time used for processing:", time_parallel[3], "seconds. \n")

simulation_data = as.data.frame(simulation_data)
colnames(simulation_data) <- c("b_tsls", "b_liml", "var_x", "var_u", "rsq", "var_v", "cov_uv", "wa_qz", "wa_delta", "tsls_se")

# 3 different approaches to normal approximation

# First, simulate using estimates from data
simu_norm_var = mean(simulation_data$var_u) / (mean(simulation_data$var_x) * mean(simulation_data$rsq))
# But, we have to shrink by dividing by sqrt(n), which divides variance by n
simu_asym_var = simu_norm_var / 100
simulation_data$normal_data = rnorm(n_sims, 0, simu_norm_var)

# Second, consider standard errors from fitting TSLS model
simulation_data$normal_med = rnorm(n_sims, 0, median(simulation_data$tsls_se))

# Finally, use N(0, 4) from analytic calculation

# 2 Approaches for weak asymptotic
# First, data. We estimate v-cov matrix and QZ*delta
sig_hat = matrix(c(mean(simulation_data$var_u), mean(simulation_data$cov_uv),
                   mean(simulation_data$cov_uv), mean(simulation_data$var_v)),
                 nrow = 2, ncol = 2)
qz_delt = mean(simulation_data$wa_qz) * mean(simulation_data$wa_delta)

# And we can simulate using the weak asymptotic distribution
uv_draws = rmvnorm(n_sims, sigma = sig_hat)

simulation_data$wa_app = uv_draws[,1] / (qz_delt + uv_draws[,2])


# Second, analytical given we know everything
true_sig = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
true_delt_qz = 0.5

exact_draws = rmvnorm(n_sims, sigma = true_sig)
simulation_data$wa_exact =  exact_draws[,1] / (0.5 + exact_draws[,2])


ggplot(simulation_data) +
  geom_density(aes(x = b_tsls, color = "Exact Distribution")) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 2), aes(color = "Normal Approximation")) +
  geom_density(aes(x = wa_app, color = "Weak Instrument Asymptotics Approximation")) +
  scale_x_continuous(name = "Distributions of 2SLS Estimator", limits = c(-6, 6)) +
  scale_y_continuous(name = "") +
  #scale_color_manual(values = c("green", "blue", "red"), labels = c("Exact Distribution", "Normal Approximation", "Weak Instrument Asymptotics Approximation")) +
  theme_light() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("out/tsls_asymptotics.png")

ggplot(simulation_data) +
  geom_density(aes(x = b_liml, color = "Exact Distribution")) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 2), aes(color = "Normal Approximation")) +
  geom_density(aes(x = wa_app, color = "Weak Instrument Asymptotics Approximation")) +
  scale_x_continuous(name = "Distributions of LIML Estimator", limits = c(-6, 6)) +
  scale_y_continuous(name = "") +
  #scale_color_manual(values = c("green", "blue", "red"), labels = c("Exact Distribution", "Normal Approximation", "Weak Instrument Asymptotics Approximation")) +
  theme_light() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("out/liml_asymptotics.png")

ggplot(simulation_data) +
  geom_density(aes(x = wa_exact, color = "Exact Distribution")) +
  geom_density(aes(x = wa_app, color = "Data Approach")) +
  scale_x_continuous(name = "Distributions of 2SLS Estimator", limits = c(-6, 6)) +
  scale_y_continuous(name = "") +
  #scale_color_manual(values = c("red", "green"), labels = c("Exact Solution", "Data Approach")) +
  theme_light() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("out/wa_approximations.png")

ggplot(simulation_data) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 2), aes(color = "Analytical")) +
  geom_density(aes(x = normal_med, color = "Median Approach")) +
  scale_x_continuous(name = "Normal Approximations", limits = c(-6, 6)) +
  scale_y_continuous(name = "") +
  #scale_color_manual(values = c("green", "blue", "red"), labels = c("Data Approach", "Analytical Approach", "Median Approach")) +
  theme_light() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("out/normal_approximations_update.png")


