rm(list = ls())
library(lhs)
library(tictoc)
source("main_tax.R")

nEval <- 240
q <- 10
nFeatures <- 1000
nTrials <- 30 
M <- 5
d <- 1
n_0 <- 24 
xmin <- matrix(rep(0, d),1,d)
xmax <- matrix(rep(1, d),1,d)
theta_0 <- matrix(c(1,0.2,1e-2),3,1)
x_grid <- matrix(seq(0,1,length=100),100,1)
k <- 4
lambda <- 0.5
param1 <- -1
param2 <- 0.39
skill <- ExpectedSkill_log(k, param1, param2)
ut <- SocialWelfare(x_grid, skill, lambda, k)
opt <- max(ut)
(xopt <- x_grid[which.max(ut)])
X <- qunif(randomLHS(n_0,d),0,1) 
Y <- QueryLaborSupply(X, k, param1, param2)
theta_0 <- matrix(rep(0.5,3),3,1)
(theta <- maximize_marginal_LLH_update(X, Y, theta_0))
#out <- BayesianOptimization(q, QueryLaborSupply, nEval, X, Y, theta, param1, param2, k, xmin, xmax, nFeatures, M, lambda, 1, TRUE)

#u_t <- SocialWelfare(x_grid, skill, lambda,k)
#u_hat_rff <- PosteriorExpectedWelfare_rff_eig(x_grid, out$X, out$Y, out$theta_opt, nFeatures, lambda)
#m_hat <- PosteriorMean(x_grid, out$X, out$Y, out$theta_opt)
#m_t <- AverageHours(x_grid, k, param1, param2)
#plot(out$X, out$Y, xlim = c(0,1), ylim = c(-0.1,1))
#lines(x_grid, m_hat, lty = 1, lw = 2)
#lines(x_grid, m_t, lty = 2, lw = 2)
#lines(x_grid, u_hat_rff, lty = 1)
#lines(x_grid, u_t, lty = 2)
#points(x_grid[which.max(u_hat_rff)], max(u_hat_rff), pch = 8) 


tic()
rep <- RepeatBayesianOptimization(1, QueryLaborSupply, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_20 <- RepeatBayesianOptimization(20, QueryLaborSupply, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_40 <- RepeatBayesianOptimization(40, QueryLaborSupply, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_60 <- RepeatBayesianOptimization(60, QueryLaborSupply, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_120 <- RepeatBayesianOptimization(120, QueryLaborSupply, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
toc()

plot(rep$reg_bar, type = "l", ylim = c(0, max(rep$reg_bar,
                                              rep_20$reg_bar,
                                              rep_40$reg_bar,
                                              rep_60$reg_bar,
                                              rep_120$reg_bar)),
  xlab = "n", ylab = "Mean Average Regret")
lines(rep_20$reg_bar, lty = 2)
lines(rep_40$reg_bar, lty = 3)
lines(rep_60$reg_bar, lty = 4)
lines(rep_120$reg_bar, lty = 5)
legend("topright", legend = c("Q=1","Q=20","Q=40","Q=60","Q=120"), lty = c(1,2,3,4,5))



# k = 2
tic()
rep_k <- RepeatBayesianOptimization(1, QueryLaborSupply, nEval, X, param1, param2, 2, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_k_20 <- RepeatBayesianOptimization(20, QueryLaborSupply, nEval, X, param1, param2, 2, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_k_40 <- RepeatBayesianOptimization(40, QueryLaborSupply, nEval, X, param1, param2, 2, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_k_60 <- RepeatBayesianOptimization(60, QueryLaborSupply, nEval, X, param1, param2, 2, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
rep_k_120 <- RepeatBayesianOptimization(120, QueryLaborSupply, nEval, X, param1, param2, 2, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, FALSE)
toc()

plot(rep_k$reg_bar, type = "l", ylim = c(0, max(rep_k$reg_bar,
                                              rep_k_20$reg_bar,
                                              rep_k_40$reg_bar,
                                              rep_k_60$reg_bar,
                                              rep_k_120$reg_bar)),
     xlab = "n", ylab = "Mean Average Regret")
lines(rep_k_20$reg_bar, lty = 2)
lines(rep_k_40$reg_bar, lty = 3)
lines(rep_k_60$reg_bar, lty = 4)
lines(rep_k_120$reg_bar, lty = 5)
legend("topright", legend = c("Q=1","Q=20","Q=40","Q=60","Q=120"), lty = c(1,2,3,4,5))











