#### Question 3
library(ggplot2)
library(ivreg)

# Read in data
ajr <- read.table("data/ajr.txt", header = T)

# OLS model
OLS_mod = lm(GDP ~ Exprop + Latitude, data = ajr)

# store estimates 
OLS_beta = summary(OLS_mod)$coefficients[2,1]
OLS_se = summary(OLS_mod)$coefficients[2,2]

OLS_est = c(OLS_beta, OLS_se)

# 2SLS model
TSLS_mod = ivreg(GDP ~ Latitude + Exprop | Mort + Latitude, data = ajr)

# store estimates
TSLS_beta = summary(TSLS_mod)$coefficients[3,1]
TSLS_se = summary(TSLS_mod)$coefficients[3,2]

TSLS_est = c(TSLS_beta, TSLS_se)

part1_out = cbind(OLS_est, TSLS_est)

# Report estimates
stargazer(part1_out, out = "out/ajr_ols_tsls.tex")

# AR statistic function
AR_calc <- function(Z, Y, X, beta) {
  n = length(Y)
  
  PZ <-  Z %*% solve(t(Z) %*% Z) %*% t(Z) 
  num <- t(Y - X * beta) %*% PZ %*% (Y - X * beta)
  denom <- t(Y - X * beta) %*% (Y - X * beta)
  AR_stat <- num/(denom/n)
  
  return(AR_stat)
}

# Setup grid
ar_grid = seq(from = 0.10, to = 2.25, by = 0.01)

# Calculate AR statistic on the grid
ar_stat = rep(NA, length(ar_grid))
for (i in 1:length(ar_grid)) {
  ar_stat[i] = AR_calc(Z = ajr$Mort, Y = ajr$GDP, X = ajr$Exprop, beta = ar_grid[i])  
}

grid_stat = as.data.frame(cbind(ar_grid, ar_stat))

# Plot AR stats
ggplot(grid_stat, aes(x = ar_grid, y = ar_stat)) + 
  geom_point() + 
  theme_light() +
  scale_x_continuous(name = "Beta") +
  scale_y_continuous(name = "AR(Beta)")
ggsave("out/ar_stats.png", width = 6, height = 3.6)

# In AR range
grid_stat$not_rej = ifelse(grid_stat$ar_stat <= qchisq(0.95, 1), 1, 0)

non_rej_b = filter(grid_stat, not_rej == 1)
ar_min = min(non_rej_b$ar_grid)
ar_max = max(non_rej_b$ar_grid)

# AR CI
ar_center = (ar_min + ar_max) / 2
ar_se = (ar_max - ar_min) / (2 * qnorm(0.975))

# Calculate confidence intervals
ar_ci = paste0("(", round(ar_min, 3), ", ", round(ar_max, 3), ")")
tsls_ci = paste0("(", round(TSLS_beta - (qnorm(0.975) * TSLS_se), 3), ", ", round(TSLS_beta + (qnorm(0.975) * TSLS_se), 3), ")")
ols_ci = paste0("(", round(OLS_beta - (qnorm(0.975) * OLS_se), 3), ", ", round(OLS_beta + (qnorm(0.975) * OLS_se), 3), ")")

cis = c(ols_ci, tsls_ci, ar_ci)
betas = c(round(OLS_beta, 3), round(TSLS_beta, 3), round(ar_center, 3))
ses = c(round(OLS_se, 3), round(TSLS_se, 3), round(ar_se, 3))

part3_out = cbind(betas, ses, cis)

# Report estimates
stargazer(part3_out, out = "out/confidence_intervals.tex")
