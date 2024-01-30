#### Question 2
library(mvtnorm)
library(ivmodel)
library(fastmatrix)

beta = 0
delta = 0.5
n = 100
z_i = rnorm(n)
sig_mat = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)

uv = rmvnorm(n, sigma = sig_mat)

x_i = (delta/(sqrt(n))) * z_i + uv[,2]
y_i = beta * x_i + uv[,1]

Z_ = as.matrix(z_i)
X_ = as.matrix(x_i)

QZ = (Z_ %*% t(Z_)) / n
QZX = (t(Z_) %*% X_) /sqrt(n)

WU_dist = kronecker.prod(sig_mat, QZ)

test = t(QZX) %*% solve(QZ)

wi_dist_draw = solve(t(QZX) %*% solve(QZ) %*% QZX) %*% t(QZX) %*% solve(QZ)

ivmodel(Y = y_i, D = x_i, Z = z_i)