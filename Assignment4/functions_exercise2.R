kf_logLik_dt <- function(par, df) {
  # par: vector of parameters
  # df: data frame with observations and inputs as columns (Y, Ta, S, I)
  # par: Could be on the form c(A11, A12, A21, A22, B11, B12, B21, B22, Q11, Q12, Q22)
  A   <- matrix() # transition matrix
  B   <- matrix() # input matrix
  Qlt <- matrix() # lower-triangle of system covariance matrix
  Q   <- Qlt %*% t(Qlt) # THAT IS!!! The system covariance matrix is given by Qlt %*% t(Qlt) (and is this symmetric positive definite)
  H   <- matrix() # observation matrix NOTE  THIS IS C ! 
  X0  <- matrix() # initial state
  R_obs <- matrix() # observation noise covariance matrix
  obs_cols <- c("Y") # observation column names
  input_cols <- c("Ta","S","I") # input column names

  # pull out data
  Y  <- as.matrix(df[, obs_cols])     # m×T
  U  <- as.matrix(df[, input_cols])   # p×T
  Tn <- nrow(df)

  # init
  n      <- nrow(A)
  x_est  <- matrix(Y[1,], n, 1)            # start state from first obs
  x_est  <- X0 
  P_est  <- diag(1e1, n)                   # X0 prior covariance
  logLik <- 0
  

  for (t in 1:Tn) {
    # prediction step
    x_pred <- # write the prediction step
    P_pred <- # write the prediction step

      # the prediction step
    if (t == 1) {
      x_pred[t] <- X0
      P_pred[t] <- P_est
    } else {
      x_pred <-  A %*% x_est + B %*% U[t, ] # NOTE x_est is x_t-1 as is was initiated last cycle 
      P_pred <-  A %*% P_est %*% t(A) + Q # 10.77 in the book P_pred is Sigmaxx,t|t
      }
  
  
    # innovation step
    y_pred <- H %*% x_pred # # NOTE C = H and then the rest is from the y= which was given
    S_t <- H %*% P_pred %*% t(H) + R_obs # (10.146) in the book 
    innovation_var <- S_t
    innov  <- Y[t, ] - y_pred # this is (10.145) from the book 

    # log-likelihood contribution
    logLik <- logLik - 0.5*(sum(log(2*pi*S_t)) + t(innov) %*% solve(S_t, innov))

    # update step
    K_t    <- P_pred * t(Q)* solve(innovation_var) # Kalman gain
    x_est  <- x_pred + K_t * (Y[t, ] - H %*% x_pred) # (10.73) from the book  reconstructed state
    P_est  <- P_pred - K_t %*% innovation_var %*% t(K_t) # (10.74) in the book reconstructed covariance
  } 

  as.numeric(logLik)
}

# Optimizer wrapper
estimate_dt <- function(start_par, df, lower=NULL, upper=NULL) {
  negLL <- function(par) -kf_logLik_dt(par, df)
  optim(
    par    = start_par, fn = negLL,
    method = "L-BFGS-B",
    lower  = lower, upper = upper,
    control= list(maxit=1000, trace=1)
  )
}

