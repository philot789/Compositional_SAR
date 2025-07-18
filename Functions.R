simulate_spatiotemporal_multivariate_p <- function(n, p, t, k, W, Psi, Pi, Beta, sig_eps = 1){
  
  Y       <- array(, dim = c(n, p, t))
  burn_in <- 100
  X       <- array(rnorm(n*p*t*k, mean = 0, sd = 1), dim = c(n, p, t, k))
  X[,,,1] <- 1 # intercept
  eps     <- array(rnorm(n*p*t, sd = sig_eps), dim = c(n, p, t))
  
  dimW <- dim(W)
  if(dimW[1] != n | dimW[2] != n){
    stop("Dimension of W is wrong")
  }
  
  vec <- function(x){
    return(as.vector(x))
  }
  
  S         <- diag(n*p) - t(Psi) %x% W
  S_inv     <- solve(S)
  
  # burn-in
  Y_burnin_t  <- array(0, dim = c(n, p))
  for(i.t in 1:burn_in){
    X_burnin_t <- array(rnorm(n*p*1*k, mean = 30, sd = 10), dim = c(n, p, 1, k))
    X_burnin_t[,,,1] <- 1
    constant <- Reduce("+", lapply(1:k, function(x) X_burnin_t[,,,x] * matrix(Beta[x,], n, p, byrow = TRUE))) # regressors
    Y_burnin_t <- array(
      S_inv %*% ( vec(constant) + 
                    vec((Y_burnin_t - W %*% Y_burnin_t %*% Psi) %*% Pi) + # AR(1) lag
                    vec(array(rnorm(n*p), dim = c(n,p))) ),  # random error
      dim = c(n, p)) # reshape to array
  }
  
  # simulations 
  
  for(i.t in 1:t){
    constant <- Reduce("+", lapply(1:k, function(x) X[,,i.t,x] * matrix(Beta[x,], n, p, byrow = TRUE))) # regressors
    if(i.t == 1){
      Y[,,i.t] <- array(
        S_inv %*% ( vec(constant) + 
                      vec((Y_burnin_t - W %*% Y_burnin_t %*% Psi) %*% Pi) + # AR(1) lag
                      vec(eps[,,i.t])),  # random error
        dim = c(n, p)) # reshape to array
    } else {
      Y[,,i.t] <- array(
        S_inv %*% ( vec(constant) + 
                      vec((Y[, , i.t - 1] - W %*% Y[, , i.t - 1] %*% Psi) %*% Pi) + # AR(1) lag
                      vec(eps[,,i.t])),  # random error
        dim = c(n, p)) # reshape to array
    }
  }
  
  ## test
  
  # res_back <- residuals(c(vec(Psi), vec(Pi), vec(Beta), 1), Y, W, X, tfe = FALSE)
  # round(eps, 4) == round(res_back, 4)
  
  return(list(Y = Y, 
              X = X))
  
}

vec <- function(x){
  return(as.vector(x))
}

qml_spatiotemporal_multivariate_p <- function(Y, W, X){
  
  dimY <- dim(Y)
  n    <- dimY[1]
  p    <- dimY[2]
  t    <- dimY[3]
  k    <- dim(X)[4] # number of regressors
  
  dimW <- dim(W)
  if(dimW[1] != n | dimW[2] != n){
    stop("Dimension of W is wrong")
  }
  
  LogLikelihood <- function(pars, Y, W, X){
    
    dimY <- dim(Y)
    n    <- dimY[1]
    p    <- dimY[2]
    t    <- dimY[3] 
    s    <- t - 1
    k    <- dim(X)[4] # number of regressors
    
    Psi     <- matrix(pars[1:(p^2)], p, p) 
    Pi      <- matrix(pars[(p^2+1):(2*p^2)], p, p) 
    Beta    <- matrix(pars[(2*p^2+1):(2*p^2 + k*p)], k, p)
    sig_u   <- pars[2*p^2 + k*p + 1]
    
    vec <- function(x){
      return(as.vector(x))
    }
    
    # start <- Sys.time()
    
    S         <- Matrix(diag(n*p) - t(Psi) %x% W)
    log_det_S <- determinant(S, logarithm = TRUE)$modulus
    
    # approximation for speed-up
    # log_det_S <- determinant(diag(n) -  Psi[1,1] * W, logarithm = TRUE)$modulus +
    #  determinant(diag(n) -  Psi[1,2] * W, logarithm = TRUE)$modulus +
    #  determinant(diag(n) -  Psi[2,1] * W, logarithm = TRUE)$modulus +
    #  determinant(diag(n) -  Psi[2,2] * W, logarithm = TRUE)$modulus

    # end <- Sys.time()
    # cat("det:", end - start, " ")
    
    # start <- Sys.time()
    
    sum_eps_2 <- 0
    for(i.t in 2:t){
      constant <- Reduce("+", lapply(1:k, function(x) X[,,i.t,x] * matrix(Beta[x,], n, p, byrow = TRUE))) # regressors
      vec_u_t  <- S %*% vec(Y[,,i.t]) - vec(constant) - vec((Y[, , i.t-1] - W %*% Y[, , i.t-1] %*% Psi) %*% Pi )
      sum_eps_2 <- sum_eps_2 + sum(vec_u_t^2)
    }
    
    # end <- Sys.time()
    # cat("sum:", end - start, " ")
    
    
    timept <- t - 1
    
    LL <- -1/2 * log(2*pi) - (timept * n * sig_u^2) / (2*p) + timept/p * log_det_S - 1/(2*p*sig_u) * sum_eps_2
    # LL <- -(timept*n*p)/2 * log(2*pi) - (n * log(sig_u^2)) / (2*p) + timept/(n*p) * log_det_S - 1/(2*n*p*sig_u^2) * sum_eps_2
    
    # cat(LL, " ")
    
    return((-1) * LL) # ((t-1)*n*p) *
  }
  
  residuals <- function(pars, Y, W, X){
    
    vec <- function(x){
      return(as.vector(x))
    }
    
    dimY <- dim(Y)
    n    <- dimY[1]
    p    <- dimY[2]
    t    <- dimY[3] 
    s    <- t - 1
    k    <- dim(X)[4] # number of regressors
    
    Psi     <- matrix(pars[1:(p^2)], p, p) 
    Pi      <- matrix(pars[(p^2+1):(2*p^2)], p, p) 
    Beta    <- matrix(pars[(2*p^2+1):(2*p^2 + k*p)], k, p)
    
    S         <- diag(n*p) - t(Psi) %x% W
    
    U_t <- array(, dim = c(n, p, t))
    for(i.t in 2:t){
      constant <- Reduce("+", lapply(1:k, function(x) X[,,i.t,x] * matrix(Beta[x,], n, p, byrow = TRUE))) # regressors
      vec_u_t  <- S %*% vec(Y[,,i.t]) - vec(constant) - vec((Y[, , i.t-1] - W %*% Y[, , i.t-1] %*% Psi) %*% Pi)
      U_t[,,i.t] <- matrix(vec_u_t, n, p)
    }
    
    return(U_t)
  }
  
  start_Psi   <- matrix(rep(0.3, p^2), p, p) 
  start_Pi    <- matrix(rep(0.3, p^2), p, p) 
  start_Beta  <- matrix(rep(0.3, k*p), k, p)
  
  start_pars <- c(vec(start_Psi), vec(start_Pi), vec(start_Beta), 1)
  
  LB_Psi <- matrix(rep(0, p^2), p, p) 
  LB_Pi  <- matrix(rep(-1, p^2), p, p)
  LB_Beta  <- matrix(rep(-1000, k*p), k, p)
  
  LB <- c(vec(LB_Psi), vec(LB_Pi), vec(LB_Beta), 0.00001)
  
  UB_Psi <- matrix(rep(1, p^2), p, p) 
  UB_Pi  <- matrix(rep(1, p^2), p, p) 
  UB_Beta  <- matrix(rep(1000, k*p), k, p)
  
  UB <- c(vec(UB_Psi), vec(UB_Pi), vec(UB_Beta), 10000)
  
  # LogLikelihood(start_pars, Y, W)
  
  out <- solnp(start_pars, fun = LogLikelihood, Y = Y, W = W, X = X, control = list(trace = TRUE),
               LB = LB, UB = UB)
  
  res <- residuals(out$pars, Y = Y, W = W, X = X)
  
  Psi_est     <- matrix(out$pars[1:(p^2)], p, p) 
  Pi_est      <- matrix(out$pars[(p^2+1):(2*p^2)], p, p) 
  Beta_est    <- matrix(out$pars[(2*p^2+1):(2*p^2 + k*p)], k, p)
  sig_u_est   <- out$pars[2*p^2 + k*p + 1]
  
  return(list(Psi_est = Psi_est, 
              Pi_est = Pi_est, 
              sig_u_est = sig_u_est, 
              Beta_est = Beta_est, 
              H = out$hessian, 
              pars_est = out$pars, 
              ll = out$values[length(out$values)], 
              residuals = res,
              n = n,
              p = p,
              t = t))
  
}

