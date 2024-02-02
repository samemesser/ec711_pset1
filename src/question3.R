#### Question 3


data_ajr = read.table('ajr.txt', header=TRUE)

#1. Estimate the effect of institutions X on log GDP per capita Y 
#   controlling for distance from the equator W by OLS
OLS_ajr <- lm(data_ajr$GDP ~ data_ajr$Exprop + data_ajr$Latitude)
#store estimate and SE 
OLS_beta <- summary(OLS_ajr)$coefficients[2,1]
OLS_se <- summary(OLS_ajr)$coefficients[2,2]

#Estimate the effect by 2SLS using settler mortality as an instrument.
twosls_ajr_1_1 <- ivreg(data_ajr$GDP ~ data_ajr$Latitude + data_ajr$Exprop | data_ajr$Mort + data_ajr$Latitude )
#store estimate, SE, and CI
TSLS_beta <- summary(twosls_ajr_1_1)$coefficients[3,1]
TSLS_se <- summary(twosls_ajr_1_1)$coefficients[3,2]
min_ci_norm <- TSLS_beta - TSLS_se*1.96
max_ci_norm <- TSLS_beta + TSLS_se*1.96

#Report estimates
stargazer(OLS_ajr, twosls_ajr_1_1)

#2. Compute the AR statistic at a grid of values between 
#   .10 and 2.25 and plot it in a figure 

#Define sequence of betas 
grid <- data.frame(cbind(matrix(seq(from = .10, to = 2.25, by = .01)),0,0))
N <- length(grid$X1)

#Define the AR statistic
AR_calc <- function(Z,Y,X,beta) {
  projection_matrix <-  Z %% solve(t(Z)%%Z)%*% t(Z) 
  numerator <- t(Y - X*beta) %% projection_matrix %% (Y - X*beta)
  denominator <- t(Y - X*beta) %*% (Y - X*beta)
  AR_stat <- numerator/(denominator/64)
  return(AR_stat)
}

#Calculate the AR stat at each value of beta
for(n in c(1:N)) {
  grid$X2[n] <- AR_calc(Z = data_ajr$Mort, Y = data_ajr$GDP, X = data_ajr$Exprop, beta=grid[n,1])
}

#plot AR stats
plot(grid$X1, grid$X2, xlab = "betas", ylab = "Anderson-Rubin Statistics")

#Find AR confidence interval 
cv <- 3.841
grid$X3 <- ifelse(grid$X2<=cv,1,0)
filtered_grid <- filter(grid, X3==1)
min_ci <- min(filtered_grid$X1)
max_ci <- max(filtered_grid$X1)

#Compute AR Estimate and SEs
AR_beta <- (max_ci + min_ci)/2
AR_se <- (max_ci - min_ci)/(2*1.96)

#3. Make tables to output 
ajr_output_CIs <- data.frame(CI_95_pct = c("Lower Bound", "Upper Bound"),
                             normal_approximation = c(min_ci_norm, max_ci_norm),
                             Anderson_Rubin = c(min_ci, max_ci))
xtable(t(ajr_output_CIs))

ajr_output_ests <- data.frame(est = c("OLS", "2SLS", "AR"),
                              estimate = c(OLS_beta,TSLS_beta,AR_beta),
                              standard_errors = c(OLS_se,TSLS_se,AR_se))
xtable(t(ajr_output_ests))