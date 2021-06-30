library(memoise)

# X1 (n x d) input matrix
# X2 (m x d) input matrix
# theta (d+2 x 1) vector of hyperparameters. Lengthscale for each dimension d
# and prior variance sigma_0^2 and noise variance sigma_n^2.
BuildCovarianceMatrix = function(X1, X2, theta){
  X1 <- as.matrix(X1) 
  X2 <- as.matrix(X2)
  n <- dim(X1)[1]
  d <- dim(X1)[2]
  m <- dim(X2)[1]
  K <- matrix(0, n, m)
  if(d==1){
    Lambda_inv <- theta[1:d]^(-2) 
  } else {
    Lambda_inv <- diag(theta[1:d]^(-2))  
  }
  for(i in 1:n){
    for(j in 1:m){
      x_diff <- matrix(X1[i,]-X2[j,],1,d)
      r2 <- x_diff %*% Lambda_inv %*% t(x_diff)
      K[i,j] <- theta[d+1]*exp(-0.5*r2)
    }
  }
  return(K)
}

# X (n x d) points sampled
# theta (d+2 x 1) vector of hyperparameters
# I have borrowed line 32-54 from Max Kasy 
posterior_prep = function(X, Y, theta){
  X <- as.matrix(X)
  C <- BuildCovarianceMatrix(X, X, theta)
  K = C + diag(theta[dim(X)[2]+2], nrow(X))
  cho = chol(K)
  Kinv <- chol2inv(cho)
  logdetK = 2*sum(log(diag(cho)))
  KinvY <- Kinv%*%Y
  list(C = C, K=K, Kinv=Kinv, KinvY=KinvY, logdetK = logdetK)
}

memo_posterior_prep = memoise(posterior_prep)

# x (q x d) points to be predicted 
# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
PosteriorMean = function(x, X, Y, theta){
  p_prep = memo_posterior_prep(X, Y, theta)
  Cs <- BuildCovarianceMatrix(x, X, theta)
  mu <- Cs%*%p_prep$KinvY
  return(mu)
}

# x (q x d) points to be predicted 
# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
# nFeatures (1 x 1) number of draws used in RFF approximation 
PosteriorMean_rff_eig = function(x, X, Y, theta, nFeatures){
  n <- dim(X)[1]
  d <- dim(X)[2]
  q <- dim(x)[1]
  X <- t(X)
  x <- t(x)
  
  W <- matrix(rnorm(nFeatures*d), nFeatures, d)*rep(theta[1:d]^(-1), nFeatures)
  b <- 2*pi*runif(nFeatures)
  
  Phi_t <- sqrt(2*theta[d+1]/nFeatures)*cos(W%*%X+matrix(rep(b,n), nFeatures, n)) # Phi_t is t(Phi) i.e. transpose of Phi 
  A <- t(Phi_t)%*%Phi_t+diag(n)*theta[d+2]
  Ainv <- chol2inv(chol(A))
  z <- Phi_t%*%Y/theta[d+2]
  m <- z - Phi_t%*%(Ainv%*%(t(Phi_t)%*%z))
  
  phi <- (sqrt(2*theta[d+1]/nFeatures)*cos(W%*%x+matrix(b, nFeatures, q)))
  mu <- t(phi)%*%m
  return(mu)
}

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
# Borrowed from Max Kasy 
marginal_llh = function(X, Y, theta){
  p_prep = memo_posterior_prep(X, Y, theta)
  ll <- -0.5*(t(Y)%*%p_prep$KinvY + p_prep$logdetK)
  return(ll)
}

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters
marginal_llh_grad = function(X, Y, theta){
  d <- dim(X)[2]
  n <- dim(X)[1]
  p_prep = memo_posterior_prep(X, Y, theta)
  C= p_prep$C
  K = p_prep$K
  Kinv = p_prep$Kinv
  KinvY = p_prep$KinvY
  
  # dK_s/dtheta_i 
  dK_sdl <- rep(list(matrix(0, n, n)),d)
  for(k in 1:d){
    for(i in 1:n){
      for(j in 1:n){
        dK_sdl[[k]][i,j] <- C[i,j]*(X[i,k]-X[j,k])^2/(theta[k]^3)
      }
    } 
  }
  dK_sdsigma2_f <- C/theta[d+1]
  dK_sdsigma2_n <- diag(n)
  
  # KinvY KinvY'-Kinv 
  KinvY_Kinv <- KinvY%*%t(KinvY)-Kinv 
  
  # dll/dtheta_i 
  grad_ll <- matrix(0, d+2, 1)
  for(k in 1:d){
    grad_ll[k] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdl[[k]]))
  }
  grad_ll[d+1] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdsigma2_f))
  grad_ll[d+2] <- 0.5*sum(diag(KinvY_Kinv%*%dK_sdsigma2_n))
  return(grad_ll)
}

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# theta (d+2 x 1) vector of hyperparameters (initial)
maximize_marginal_LLH_update = function(X, Y, theta){
  d <- dim(X)[2]
  theta_min <- rep(1e-3, d+2) 
  theta_max <- rbind(matrix(rep(1e2, d),d,1), 
                     matrix(var(Y),2,1))
  
  res <- optim(par = theta,
               fn = function(theta) -marginal_llh(X, Y, theta),
               gr = function(theta) -marginal_llh_grad(X, Y, theta), 
               lower = theta_min,
               upper = theta_max,
               method = "L-BFGS-B",
               control = list(maxit = 10, trace=0)) 
  
  theta_opt <- res$par
  return(theta_opt)
} 

# X (n x d) points sampled
# Y (n x 1) value at points sampled
# M (1 x 1) number of restarts 
maximize_marginal_LLH = function(X, Y, M){
  d <- dim(X)[2]
  obj_list = matrix(0, 1, M)
  theta_list = matrix(0, d+2, M)
  theta_min <- rep(1e-3, d+2) # was previously 1e-6
  theta_max <- rbind(matrix(rep(10, d),d,1), matrix(rep(var(Y), 2),2,1)) # was 1 and 10
  theta_0 <- matrix(runif((d+2)*M, theta_min[1], theta_max[1]), d+2, M) # try restricting this 
  theta_0[2,] <- rep(var(Y), M) # note 
  
  for(i in 1:M){
    res <- optim(par = theta_0[,i],
                 fn = function(theta) -marginal_llh(X, Y, theta),
                 gr = function(theta) -marginal_llh_grad(X, Y, theta), 
                 lower = theta_min,
                 upper = theta_max,
                 method = "L-BFGS-B",
                 control = list(trace = 0))
    theta_opt <- res$par
    ll_opt <- res$value 
    theta_list[,i] <- theta_opt 
    obj_list[i] <- ll_opt 
  }
  ind <- which.min(obj_list)
  best <- theta_list[, ind]
  best <- matrix(best,d+2,1)
  list(theta_list = theta_list,
       obj_list = obj_list,
       best = best 
  )
} 

ThompsonSpectralSampling = function(X, Y, theta, xmin, xmax, nFeatures, M, lambda){
  n <- dim(X)[1]
  d <- dim(X)[2]
  X <- t(X)
  xmin <- t(xmin)
  xmax <- t(xmax)
  
  W <- matrix(rnorm(nFeatures*d), nFeatures, d)*rep(theta[1:d]^(-1), nFeatures)
  b <- 2*pi*runif(nFeatures)
  noise <- rnorm(nFeatures)
  
  Phi_t <- sqrt(2*theta[d+1]/nFeatures)*cos(W%*%X+matrix(rep(b,n), nFeatures, n)) # Phi_t is t(Phi) i.e. transpose of Phi 
  A <- t(Phi_t)%*%Phi_t+diag(n)*theta[d+2]
  Ainv <- chol2inv(chol(A))
  z <- Phi_t%*%Y/theta[d+2]
  m <- z - Phi_t%*%(Ainv%*%(t(Phi_t)%*%z))
  eig <- eigen(A)
  R <- (sqrt(eig$values)*(sqrt(eig$values)+sqrt(theta[d+2])))^(-1)
  omega <- noise - (Phi_t%*%(eig$vectors%*%(R*(t(eig$vectors)%*%(t(Phi_t)%*%noise)))))+m
  
  obj = function(x){
    phi <- (sqrt(2*theta[d+1]/nFeatures)*cos(W%*%x+matrix(b, nFeatures, 1)))
    term1 <- sqrt(2*theta[d+1]/nFeatures)*(-lambda*(sin(W%*%x+matrix(b,nFeatures,1))-sin(matrix(b,nFeatures,1))))
    term2 <- matrix(rep(x, each = nFeatures), nFeatures, 1)*phi
    phi_u <- term1*matrix(rep(W^-1),nFeatures, 1)+term2
    f <- -t(omega)%*%phi_u
    return(f)
  }
  gradobj = function(x){
    term1 <- -lambda*(matrix(rep(cos(W%*%x+b),d), nFeatures, d)*W)
    term2 <- matrix(rep(cos(W%*%x+b),d), nFeatures, d)
    term3 <- (matrix(rep(sin(W%*%x+b),d), nFeatures, d)*W)*matrix(rep(x, d),nFeatures, d)
    nablaf <- t(omega)%*%(term1+term2+term3)*(sqrt(2*theta[d+1]/nFeatures))
    return(nablaf)
  }
  
  x_0 <- matrix(runif(d*M, xmin[1], xmax[1]),d,M)
  x_list <- matrix(0, M, d)
  y_list <- matrix(0, M, 1)
  for(i in 1:(M)){
    res <- optim(par = matrix(x_0[,i], d, 1),
                 fn = obj,
                 gr = gradobj,
                 lower = xmin,
                 upper = xmax,
                 method = "L-BFGS-B")
    x_list[i,] <- res$par
    y_list[i] <- -res$value
  }
  x_next <- matrix(x_list[which.max(y_list),],1,d)
  y_best <- max(y_list)
  list(x_next = x_next,
       y_best = y_best)
}

MaximizePosteriorExpectedWelfare = function(X, Y, theta, xmin, xmax, M, nFeatures, lambda){
  obj_list <- matrix(0, M, 1)
  x_list <- matrix(0, M, d)
  x_0 <- matrix(runif(d*M, xmin, xmax), M, d)
  for(i in 1:M){
    res <- optim(par = x_0[i,], 
                 fn = function(x) -PosteriorExpectedWelfare_rff_eig(t(x), X, Y, theta, nFeatures, lambda),
                 method = "L-BFGS-B",
                 upper = xmax,
                 lower = xmin, 
                 control = list(trace = 0))
    obj_list[i] <- res$value
    x_list[i,] <- res$par
  }
  x_best <- x_list[which.min(obj_list),]
  y_best <- -min(obj_list)
  list(x_best = x_best,
       y_best = y_best)
} 

ExpectedSkill_log = function(k, mu, sigma){
  integrand <- function(x) {x^(1/k)*(1/(x*sigma*sqrt(2*pi)))*exp(-(log(x)-mu)^2/(2*sigma^2))}
  out <- integrate(integrand, lower = 0, upper = Inf)
  return(out$value)
}

SocialWelfare = function(t, expectedskill, lambda = 0.5, k = 2){
  ut <- (-lambda*(k/(1+k))*(1-(1-t)^((1+k)/k))+t*(1-t)^(1/k))*expectedskill 
  return(ut)
}

QueryLaborSupply = function(t, k = 2, param1 = 0, param2 = 1){
  n <- rlnorm(length(t), param1, param2) 
  y <- (n*(1-t))^(1/k)
  return(y)
}

AverageHours = function(x, k = 2, param1 = 0, param2 = 1) {
  skill <- ExpectedSkill_log(k, param1, param2)
  m = skill*(1-x)^(1/k)
  return(m)
}

BayesianOptimization = function(q, FUN, nEval, X, Y, theta, param1, param2, k, xmin, xmax, nFeatures, M, lambda, trial = 1, plt = FALSE){
  theta_opt <- theta
  for(i in 1:(as.integer(nEval/q))){
    x_batch <- matrix(0, q, d)
    y_batch <- matrix(0, q, 1)
    for(j in 1:q){
      ss <- ThompsonSpectralSampling(X, Y, theta_opt, xmin, xmax, nFeatures, M, lambda)
      y_batch[j] <- as.matrix(FUN(as.matrix(ss$x_next), k, param1, param2))
      x_batch[j] <- ss$x_next
      cat(paste("Trial:", trial, "\n"))
      cat(paste("Iteration:", j, "\n"))
      cat(paste("Batch:", i, "\n"))
      cat(paste("Y:", y_batch[j], "\n"))
      cat(paste("X", "current:", x_batch[j],"\n"))
      cat(paste("l", "current:", theta_opt[1],"\n"))
      cat(paste("sigma_0^2", "current:", theta_opt[2],"\n"))
      cat(paste("sigma_n^2", "current:", theta_opt[3],"\n"))
      if(plt == TRUE){
        u_t <- SocialWelfare(x_grid, skill, lambda, k)
        u_hat_rff <- PosteriorExpectedWelfare_rff_eig(x_grid, X, Y, theta_opt, nFeatures, lambda)
        m_hat <- PosteriorMean(x_grid, X, Y, theta_opt)
        m_t <- AverageHours(x_grid, k, param1, param2)
        plot(X, Y, xlim = c(0,1), ylim = c(-0.1,1))
        lines(x_grid, m_hat, lty = 2, lw = 2)
        lines(x_grid, m_t)
        lines(x_grid, u_hat_rff, lty = 2, lw = 2, col = "blue")
        lines(x_grid, u_t, col = "blue")
        points(ss$x_next, ss$y_best, col = "red", pch = 8)
        points(x_grid[which.max(u_hat_rff)], max(u_hat_rff), col = "green", pch = 8) 
      }   
    }
    X <- rbind(X, x_batch) 
    Y <- rbind(Y, y_batch)
    if(q==1){
      div = 10
    } else {
      div = 1
    }
    if(i%%div==0){ 
      if(i!=(nEval/q)){
        theta_opt <- maximize_marginal_LLH_update(X, Y, theta_opt)
      }
    }
  }
  Y_u <- PosteriorExpectedWelfare_rff_eig(X, X, Y, theta_opt, nFeatures, lambda)
  x_opt <- MaximizePosteriorMean(X, Y, theta, xmin, xmax, M, nFeatures)
  list(X = X,
       Y = Y,
       Y_u = Y_u,
       theta_opt = theta_opt,
       x_opt = x_opt)
} 

RepeatBayesianOptimization = function(q, FUN, nEval, X, param1, param2, k, xmin, xmax, nFeatures, M, lambda, nTrials, n_0, plt = FALSE){
  skill <- ExpectedSkill_log(k, param1, param2)
  ut <- SocialWelfare(x_grid, skill, lambda, k)
  opt <- max(ut)
  reg_list <- matrix(0, nEval, nTrials)
  Y <- QueryLaborSupply(X, k, param1, param2)
  theta_0 <- matrix(rep(0.5,3),3,1)
  theta <- maximize_marginal_LLH_update(X, Y, theta_0)
  for(i in 1:nTrials){
    X_temp <- BayesianOptimization(q, FUN, nEval, X, Y, theta, param1, param2, k, xmin, xmax, nFeatures, M, lambda, i, plt)$X
    reg_list[,i] <- PrepareAverageRegret(X_temp, opt, nEval, lambda, k, param1, param2, nFeatures)
  }
  reg_bar <- apply(reg_list, 1, mean)
  list(reg_list = reg_list, 
       reg_bar = reg_bar)
}

MaximizePosteriorMean = function(X, Y, theta, xmin, xmax, M, nFeatures){
  d <- dim(X)[2]
  obj_list <- matrix(0, M, 1)
  x_list <- matrix(0, M, d)
  x_0 <- matrix(runif(d*M, xmin, xmax), M, d)
  for(i in 1:M){
    res <- optim(par = x_0[i,], 
                 fn = function(x) - PosteriorExpectedWelfare_rff_eig(t(x), X, Y, theta, nFeatures),
                 method = "L-BFGS-B",
                 upper = xmax,
                 lower = xmin, 
                 control = list(trace = 0))
    obj_list[i] <- - res$value
    x_list[i,] <- res$par
  }
  x_best <- x_list[which.max(obj_list),]
  y_best <- max(obj_list)
  list(x_best = x_best,
       y_best = y_best)
} 

PosteriorExpectedWelfare_rff_eig = function(x, X, Y, theta, nFeatures, lambda = 0.5){
  n <- dim(X)[1]
  d <- dim(X)[2]
  q <- dim(x)[1]
  X <- t(X)
  x <- t(x)
  
  W <- matrix(rnorm(nFeatures*d), nFeatures, d)*rep(theta[1:d]^(-1), nFeatures)
  b <- 2*pi*runif(nFeatures)
  
  Phi_t <- sqrt(2*theta[d+1]/nFeatures)*cos(W%*%X+matrix(rep(b,n), nFeatures, n)) # Phi_t is t(Phi) i.e. transpose of Phi 
  A <- t(Phi_t)%*%Phi_t+diag(n)*theta[d+2]
  Ainv <- chol2inv(chol(A))
  z <- Phi_t%*%Y/theta[d+2]
  m <- z - Phi_t%*%(Ainv%*%(t(Phi_t)%*%z))
  
  phi <- (sqrt(2*theta[d+1]/nFeatures)*cos(W%*%x+matrix(b, nFeatures, q)))
  term1 <- sqrt(2*theta[d+1]/nFeatures)*(-lambda*(sin(W%*%x+matrix(b,nFeatures,q))-sin(matrix(b,nFeatures,q))))
  term2 <- matrix(rep(x, each = nFeatures), nFeatures, q)*phi
  phi_u <- term1*matrix(rep(W^-1),nFeatures, q)+term2
  mu_hat <- t(phi_u)%*%m
  return(mu_hat)
}

PosteriorExpectedWelfare = function(x, X, Y, theta, lambda = 0.5){
  p_prep = memo_posterior_prep(X, Y, theta)
  Kinv <- p_prep$Kinv 
  D <- matrix(0, length(x), length(X))
  for(i in 1:length(x)){
    for(j in 1:length(X)){
      D[i,j] = (theta[2]/dnorm(0))*(-lambda*theta[1]*(pnorm((x[i]-X[j])/theta[1])-pnorm((-X[j])/theta[1]))+x[i]*dnorm((x[i]-X[j])/theta[1]))
    }
  }
  u_hat <- D%*%Kinv%*%Y
  return(u_hat)
}

PlotSampleAverageRegret = function(Y, opt, nEval){
  avg_regret <- rep(0, nEval)
  for(i in 1:nEval){
    avg_regret[i] <- opt-mean(Y[1:(i+(length(Y)-nEval))]) 
  }
  plot(avg_regret, type = "l", xlab = "n", ylab = "Average regret")
}

PrepareSampleAverageRegret = function(Y, opt, nEval){
  avg_regret <- rep(0, nEval+1)
  index <- rep(0, nEval+1)
  for(i in 1:(nEval+1)){
    avg_regret[i] <- opt-mean(Y[1:(i+(length(Y)-nEval)-1)]) 
    index[i] <- i-1 
  }
  list(avg_regret = avg_regret,
       index = index)
} 

PlotAverageRegret = function(X, opt, nEval,lambda, k, param1, param2, nFeatures){
  skill <- ExpectedSkill_log(k, param1, param2)
  Y <- SocialWelfare(X, skill, lambda, k)
  avg_regret <- rep(0, nEval)
  for(i in 1:nEval){
    avg_regret[i] <- opt-mean(Y[1:(i+(length(Y)-nEval))]) 
  }
  plot(avg_regret, type = "l", xlab = "n", ylab = "Average regret")
}

PrepareAverageRegret = function(X, opt, nEval, lambda, k, param1, param2, nFeatures){
  skill <- ExpectedSkill_log(k, param1, param2)
  Y <- SocialWelfare(X, skill, lambda, k)
  avg_regret <- rep(0, nEval)
  for(i in 1:nEval){
    avg_regret[i] <- opt-mean(Y[1:(i+(length(Y)-nEval))]) 
  }
  return(avg_regret)
}




