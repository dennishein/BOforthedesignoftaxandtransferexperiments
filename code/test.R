rm(list = ls())
library(lhs)
library(tictoc)
library(latex2exp)
source("main.R") 

# Load data 
X_data <- read.csv("/Users/dennishein/Documents/BayesOpt thesis/Oxford thesis/Code/Coinsurance/X.csv")
Y_data <- read.csv("/Users/dennishein/Documents/BayesOpt thesis/Oxford thesis/Code/Coinsurance/Y.csv")
XX <- as.matrix(1-X_data[1:nrow(X_data),1])
YY <- as.matrix(Y_data[1:nrow(Y_data),1]) 

# Fit ground truth and obtain residuals 
nFeatures <- 1000
theta_data <- matrix(rep(1,3),1,3)
YYhat <- PosteriorMean_rff(XX, XX, YY, theta_data, nFeatures)
res <- YY-YYhat

# Plot mhat and uhat 
x_grid <- matrix(seq(0,1,length=100),100,1)
mhat <- PosteriorMean_rff(x_grid, XX, YY, theta_data, nFeatures)
uhat <- PosteriorExpectedWelfare_rff(x_grid, XX, YY, theta_data, nFeatures)
opt <- max(uhat) 
(topt <- x_grid[which.max(uhat)]) # optimal x
plot(x_grid, mhat, type = "l", ylim = c(0,2100),  xlab = "t", ylab = "") 
lines(x_grid, uhat, lty = 2)
points(topt, opt, pch = 8)
legend("topleft", legend=c(TeX("$\\mu_{\\tilde{N}}(t)$"), TeX("$\\upsilon_{\\tilde{N}}(t)$")),lty=1:2)

nEval <- 500
nTrials <- 50 # was 50 
q <- 1
M <- 5
d <- 1
n_0 <-  50
xmin <- matrix(rep(0, d),1,d)
xmax <- matrix(rep(1, d),1,d)
lambda <- 1.5
opt <- max(uhat)
(topt <- x_grid[which.max(uhat)])
X <- matrix(randomLHS(n_0, 1), n_0, 1) # initial sample 
version = "wild"
Y <- query_m(X, XX, YY, theta_data, nFeatures, res, version)
theta <- matrix(c(0.2,1,1), 1, 3) # 0.2,1,1


tic()
rep <- RepeatBatchBayesianOptimization(query_m, nEval, X, theta, XX, YY, theta_data, res, q, version, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, opt)
plot(rep$reg_bar, type = "l")
rep_50 <- RepeatBatchBayesianOptimization(query_m, nEval, X, theta, XX, YY, theta_data, res, 25, version, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, opt)
plot(rep_50$reg_bar, type = "l")
rep_100 <- RepeatBatchBayesianOptimization(query_m, nEval, X, theta, XX, YY, theta_data, res, 50, version, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, opt)
plot(rep_100$reg_bar, type = "l")
rep_200 <- RepeatBatchBayesianOptimization(query_m, nEval, X, theta, XX, YY, theta_data, res, 100, version, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, opt)
plot(rep_200$reg_bar, type = "l")
toc()
plot(rep$reg_bar, type = "l", ylim = c(min(rep$reg_bar, rep_50$reg_bar,  rep_100$reg_bar, rep_200$reg_bar),
                                        max(rep$reg_bar, rep_50$reg_bar, rep_100$reg_bar, rep_200$reg_bar)),
     xlab = "n", ylab = "Mean Average Regret")
lines(rep_50$reg_bar, type = "l", lty = 2)
lines(rep_100$reg_bar, lty = 3)
lines(rep_200$reg_bar, lty = 4)
legend("topright", legend = c("Q=1","Q=25","Q=50", "Q=100"), lty = c(1,2,3,4))

