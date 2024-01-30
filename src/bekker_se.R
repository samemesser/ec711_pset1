##' bekker_se.R
##'
##' This procedure calculates standard errors as in Bekker (1994)
##'
##' @title Bekker standard error calculation
##' 
##' @param ivmodel A fit IV model using the ivmodel package
##' @param Y The outcome variable
##' @param D The endogenous variable description 
##' @param Z The instruments
##' @param X The control variables (must include intercept) description
##' 
##' @return The standard error robust to many instruments
##' 
##' @author Samuel Messer
##' @export
##' 


bekker_se <- function(ivmodel, Y, D, Z, X)
{
  # Estimates from LIML
  b_tilde = as.vector(c(coefOther(ivmodel)$LIML[, "Estimate"], coef(ivmodel)["LIML","Estimate"]))
  # Vars used to estimate Y
  X_ = as.matrix(cbind(X, D))
  # Vars used to estimate X
  Z_ = as.matrix(cbind(X, Z))
  
  # Estimated errors
  u_tilde = Y - X_ %*% b_tilde
  
  # Useful quantities
  UU = as.numeric(t(u_tilde) %*% u_tilde)
  
  # Components of Bekker's SE
  sig_u_sq = UU/nrow(u_tilde)
  a_tilde = as.numeric(t(u_tilde) %*% Z_ %*% solve(t(Z_) %*% Z_) %*% t(Z_) %*% u_tilde / UU)
  X_tilde = X_ - (u_tilde %*% (t(u_tilde) %*% X_)) / UU
  H_hat = t(X_) %*% Z_ %*% solve(t(Z_) %*% Z_) %*% t(Z_) %*% X_ - a_tilde * (t(X_) %*% X_)
  Sig_hat = sig_u_sq * ((1 - a_tilde)^2 * t(X_tilde) %*% Z_ %*% solve(t(Z_) %*% Z_) %*% t(Z_) %*% X_tilde + a_tilde^2 * t(X_tilde) %*% X_tilde - a_tilde^2 * t(X_tilde) %*% Z_ %*% solve(t(Z_) %*% Z_) %*% t(Z_) %*% X_tilde)
  
  bekker_var = solve(H_hat) %*% Sig_hat %*% solve(H_hat)
  
  bekker_se = as.numeric(sqrt(diag(bekker_var))["D"])
  
  return(bekker_se)
}