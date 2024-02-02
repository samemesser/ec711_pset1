library(tidyverse)
library(readstata13)
library(fastDummies)
library(ivmodel)
library(doParallel)
library(foreach)
library(doRNG)
library(mvtnorm)
library(data.table)
library(stargazer)

rm(list = ls())

time = proc.time()

source("src/bekker_se.R")

hhn <- read.dta13("data/hhn.dta")

hhn_factors <- hhn %>%
  mutate(QOB = as.factor(QOB), STATE = as.factor(STATE), YOB = as.factor(YOB))


# Suppose instead we take last of every group to be excluded
hhn_factors$QOB <- relevel(hhn_factors$QOB, ref = tail(levels(hhn_factors$QOB), 1))
hhn_factors$YOB <- relevel(hhn_factors$YOB, ref = tail(levels(hhn_factors$YOB), 1))
hhn_factors$STATE <- relevel(hhn_factors$STATE, ref = tail(levels(hhn_factors$STATE), 1))

# This generates all of the interaction columns
interacts <- model.matrix(~ QOB + QOB:YOB + QOB:STATE, data = hhn_factors)

hhn_full <- cbind(hhn_factors, interacts)

# QOB1 is the excluded group, so we want it removed everywhere
cols_to_remove <- c(grep("QOB4", colnames(hhn_full)), 4, 7, 8)
hhn_all_inst <- hhn_full[,-cols_to_remove]
hhn_clean <- hhn_all_inst[,c(1:2, 6, 3:5, 7:186)]

# Write subset of the file so we can skip time consuming steps when we parallelize
hhn_export <- hhn_clean[,c(3:186)]
fwrite(hhn_export, file = "intermediate_data/hhn_export.csv")

# Remove extraneous datasets
rm(hhn, hhn_all_inst, hhn_factors, hhn_full, interacts, hhn_export, cols_to_remove)

#### Do some sort of validation exercise w/ Jimmy/Erin
# subse <- hhn_clean %>%
#   filter(`QOB2:YOB31` == 1) %>%
#   select(QOB, YOB, QOB2, `QOB2:YOB31`)

# The same across specifications
wage <- hhn_clean$LWKLYWGE
educ <- hhn_clean$EDUC
controls <- hhn_clean[,3:6] # Including intercept

# Different instrument specifications
instr_s1 <- hhn_clean[,7:9] # QOB only
instr_s2 <- hhn_clean[,7:36] # QOB + QOB*YOB
instr_s3 <- hhn_clean[,7:186] # QOB + QOB*YOB + QOB*STATE

# IVmodels
mod_s1 <- ivmodel(Y=wage, D=educ, Z=instr_s1, X=controls, intercept = FALSE)
mod_s2 <- ivmodel(Y=wage, D=educ, Z=instr_s2, X=controls, intercept = FALSE)
mod_s3 <- ivmodel(Y=wage, D=educ, Z=instr_s3, X=controls, intercept = FALSE)

tsls_est = c(coef(mod_s1)["TSLS","Estimate"], coef(mod_s2)["TSLS","Estimate"], coef(mod_s3)["TSLS","Estimate"])
tsls_se = c(coef(mod_s1)["TSLS","Std. Error"], coef(mod_s2)["TSLS","Std. Error"], coef(mod_s3)["TSLS","Std. Error"])
liml_est = c(coef(mod_s1)["LIML","Estimate"], coef(mod_s2)["LIML","Estimate"], coef(mod_s3)["LIML","Estimate"])
liml_se = c(coef(mod_s1)["LIML","Std. Error"], coef(mod_s2)["LIML","Std. Error"], coef(mod_s3)["LIML","Std. Error"])
bekker = c(bekker_se(mod_s1, wage, educ, instr_s1, controls), bekker_se(mod_s2, wage, educ, instr_s2, controls), bekker_se(mod_s3, wage, educ, instr_s3, controls))

base_results = rbind(tsls_est, tsls_se, liml_est, liml_se, bekker)
colnames(base_results) = c("S1", "S2", "S3")

stargazer(base_results, out = "out/results_last_excluded.tex")


######### Part (b)
# Collect parameters for simulation

# Betas
non_edg_coefs = coefOther(mod_s1)$LIML[, "Estimate"] 
b0 = non_edg_coefs[1] #Intercept
b2 = non_edg_coefs[2:4] #Control estimates
b1 = coef(mod_s1)["LIML","Estimate"] #Endogenous beta
b_ = as.matrix(c(b0, b2, b1))

# To get the parameters for the equation defining X, we run OLS with X on
# instruments and controls
# Pis
mod_x <- lm(educ ~ hhn_clean$MARRIED + hhn_clean$RACE + hhn_clean$SMSA
            + hhn_clean$QOB1 + hhn_clean$QOB2 + hhn_clean$QOB3)
p_ = as.matrix(mod_x$coefficients)

# Errors
# U
X_ = cbind(1, hhn_clean$MARRIED, hhn_clean$RACE, hhn_clean$SMSA, educ)
u_ = wage - X_ %*% b_
sig_u = sd(u_)

#V
Z_ = cbind(1, hhn_clean$MARRIED, hhn_clean$RACE, hhn_clean$SMSA,
           hhn_clean$QOB1, hhn_clean$QOB2, hhn_clean$QOB3)
v_ = educ - Z_ %*%  p_
sig_v = sd(v_)

# Cov(u, v)
sig_uv = cor(u_, v_)

# Sigma for errors
sigma_hhn = c(sig_u^2, sig_uv, sig_uv, sig_v^2)

fwrite(list(b_), file = "intermediate_data/betas.txt")
fwrite(list(p_), file = "intermediate_data/pis.txt")
fwrite(list(sigma_hhn), file = "intermediate_data/sigma_hhn.txt")

num_datasets = 1000

ncore <- detectCores()
cl <- makeCluster(10, type = "PSOCK") #ncore - 1 usually
registerDoParallel(cl)

# From random.org
set.seed(191116266)
time_parallel <- system.time({
  liml_bekker <- foreach(i = 1:num_datasets
             , .combine = 'rbind'
             , .packages = c("ivmodel", "mvtnorm", "data.table")
             ) %dorng% {

               cat("Currently processing dataset", i, "\n")

               source("src/bekker_se.R")
               hhn_export <- fread("intermediate_data/hhn_export.csv")
               betas <- as.matrix(fread("intermediate_data/betas.txt"))
               pis <- as.matrix(fread("intermediate_data/pis.txt"))
               sigma_hhn <- matrix(as.matrix(fread("intermediate_data/sigma_hhn.txt")), nrow = 2, ncol = 2)

               data_length = 329509
               errors = rmvnorm(data_length, sigma = sigma_hhn)

               Z_ = as.matrix(hhn_export[,1:7])

               X_fit = Z_ %*% pis + errors[,2]

               X_ = as.matrix(cbind(hhn_export[,1:4], X_fit))

               Y_fit = X_ %*% betas + errors[,1]

               # Name columns
               colnames(X_fit)[1] = "D"
               colnames(Y_fit)[1] = "Y"

               controls <- hhn_export[,1:4]

               # Different instrument specifications
               instr_s1 <- hhn_export[,5:7] # QOB only
               instr_s2 <- hhn_export[,5:34] # QOB + QOB*YOB
               instr_s3 <- hhn_export[,5:184] # QOB + QOB*YOB + QOB*STATE

               mod_s1 <- ivmodel(Y=Y_fit, D=X_fit, Z=instr_s1, X=controls, intercept = FALSE)
               mod_s2 <- ivmodel(Y=Y_fit, D=X_fit, Z=instr_s2, X=controls, intercept = FALSE)
               mod_s3 <- ivmodel(Y=Y_fit, D=X_fit, Z=instr_s3, X=controls, intercept = FALSE)
        
               tsls_est_s1 = coef(mod_s1)["TSLS","Estimate"]
               tsls_est_s2 = coef(mod_s2)["TSLS","Estimate"]
               tsls_est_s3 = coef(mod_s3)["TSLS","Estimate"]
               tsls_se_s1 = coef(mod_s1)["TSLS","Std. Error"]
               tsls_se_s2 = coef(mod_s2)["TSLS","Std. Error"]
               tsls_se_s3 = coef(mod_s3)["TSLS","Std. Error"]
               liml_est_s1 = coef(mod_s1)["LIML","Estimate"]
               liml_est_s2 = coef(mod_s2)["LIML","Estimate"]
               liml_est_s3 = coef(mod_s3)["LIML","Estimate"]
               liml_se_s1 = coef(mod_s1)["LIML","Std. Error"]
               liml_se_s2 = coef(mod_s2)["LIML","Std. Error"]
               liml_se_s3 = coef(mod_s3)["LIML","Std. Error"]
               bekker_s1 = bekker_se(mod_s1, Y_fit, X_fit, instr_s1, controls)
               bekker_s2 = bekker_se(mod_s2, Y_fit, X_fit, instr_s2, controls)
               bekker_s3 = bekker_se(mod_s3, Y_fit, X_fit, instr_s3, controls)
               
               out = c(tsls_est_s1, tsls_se_s1, liml_est_s1, liml_se_s1, bekker_s1, tsls_est_s2, tsls_se_s2, liml_est_s2, liml_se_s3, bekker_s2, tsls_est_s3, tsls_se_s3, liml_est_s3, liml_se_s3, bekker_s3)

               out
  }
})
stopCluster(cl)

cat("Time used for processing:", time_parallel[3], "seconds. \n")

simulated_estimates <- as.data.frame(liml_bekker)
colnames(simulated_estimates) <- c("tsls_est_s1", "tsls_se_s1", "liml_est_s1", "liml_se_s1", "bekker_s1", "tsls_est_s2", "tsls_se_s2", "liml_est_s2", "liml_se_s3", "bekker_s2", "tsls_est_s3", "tsls_se_s3", "liml_est_s3", "liml_se_s3", "bekker_s3")

write.csv(simulated_estimates, file = "intermediate_data/sim_ests.csv")

cat("All done! :) \n", "Total analysis time:", (proc.time() - time)[3], " seconds. \n")

#####
# Assume all prior was run on cluster, need to load data here!
sim_est <- read.csv("intermediate_data/sim_ests.csv")
true_beta <- b1

sim_errors <- sim_est[,c(1, 3, 6, 8, 11, 13)] - true_beta

squared_errors <- sim_errors^2

ci_coverage <- between(0, sim_errors[,c(1:6, 2, 4, 6)] - qnorm(0.975)*sim_est[,c(2, 4, 7, 9, 12, 14, 5, 10, 15)],
                          sim_errors[,c(1:6, 2, 4, 6)] + qnorm(0.975)*sim_est[,c(2, 4, 7, 9, 12, 14, 5, 10, 15)]
)

bias <- c(mean(sim_errors[,1]), mean(sim_errors[,2]), NA, mean(sim_errors[,3]),
          mean(sim_errors[,4]), NA, mean(sim_errors[,5]), mean(sim_errors[,6]), NA)
rmse <- c(sqrt(mean(squared_errors[,1])), sqrt(mean(squared_errors[,2])), NA,
          sqrt(mean(squared_errors[,3])), sqrt(mean(squared_errors[,4])), NA,
          sqrt(mean(squared_errors[,5])), sqrt(mean(squared_errors[,6])), NA)
cove <- c(mean(ci_coverage[,1]), mean(ci_coverage[,2]), mean(ci_coverage[,7]),
          mean(ci_coverage[,3]), mean(ci_coverage[,4]), mean(ci_coverage[,8]),
          mean(ci_coverage[,5]), mean(ci_coverage[,6]), mean(ci_coverage[,9]))

results_out <- cbind(bias, rmse, cove)
rownames(results_out) <-  c("tsls_s1", "liml_s1", "bekker_s1",
                            "tsls_s2", "liml_s2", "bekker_s2",
                            "tsls_s3", "liml_s3", "bekker_s3")

stargazer(results_out, out = "out/simulation_est.tex")
write.csv(results_out, "out/simulation_results_qob4.csv")

