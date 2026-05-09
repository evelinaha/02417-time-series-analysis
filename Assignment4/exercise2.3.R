kf_logLik_dt <- function(par, df) {
  # # par: vector of parameters
  # df: data frame with observations and inputs as columns (Y, Ta, S, I)
  # par: Could be on the form c(A11, A12, A21, A22, B11, B12, B21, B22, Q11, Q12, Q22)
  # A is a 2x2 matrix
  A <- matrix(par[1:4], nrow = 2, byrow = TRUE)
  # B is 2x3 (maps 3 inputs to 2 states)
  B <- matrix(par[5:10], nrow = 2, byrow = TRUE)
  # H is 1x2 (maps 2 states to 1 observation Y)
  H <- matrix(c(1, 0), nrow = 1, ncol = 2)    # C
  #H <- matrix(par[11:12], nrow = 1, ncol = 2) # C
  # Initial values for X are the first observation values (initial guess (10.79))
  X0 <- matrix(c(df$Y[1], 0), nrow = 2)
  
  obs_cols <- c("Y")
  input_cols <- c("Ta","S","I") 
  
  # System Noise Q (using lower triangular matrix Qlt to keep it positive definite)
  Qlt <- matrix(0, 2, 2)
  Qlt[1,1] <- par[11]
  Qlt[2,1] <- par[12]
  Qlt[2,2] <- par[13]
  Q <- Qlt %*% t(Qlt) # Sigma_1 
  
  # Observation Noise (scalar) and Initial State
  R_obs <- matrix(exp(par[14])) # exp ensures variance is always positive (Sigma_2)
  
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
  
  # Kalman loop
  for (t in 1:Tn) {
    # Prediction step
    # u_t is Ta, S, I at time t
    u_t <- matrix(U[t, ], 3, 1)
    
    if (t == 1) {
      x_pred <- X0
      P_pred <- diag(10, 2) # Start with some initial uncertainty
    } else {
      # x_{t+1|t} = A * x_{t|t} + B * u_t
      x_pred <- A %*% x_est + B %*% u_t
      P_pred <- A %*% P_est %*% t(A) + Q    # Sigma_xx
    }
    
    # Innovation step
    y_pred <- H %*% x_pred
    S_t    <- as.numeric(H %*% P_pred %*% t(H) + R_obs)  # innovation variance 
    innov  <- Y[t] - y_pred
    
    # Log-likelihood
    # Note: Use det(S_t) for the log-likelihood calculation
    logLik <- logLik - 0.5 * (log(S_t) + innov^2 / S_t)  # look at given code
    
    # Update step (reconstruction)
    #K_t   <- P_pred %*% t(H) %*% solve(S_t)   # Kalman gain
    K_t <- (P_pred %*% t(H)) / S_t              # Kalman gain
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(2) - K_t %*% H) %*% P_pred
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


# Read the data
setwd("C:\\Users\\Anja\\OneDrive\\Skrivebord\\Time Series Analysis\\Assignment4\\assignment4_2026")
transformer_data <- read.csv("transformer_data.csv")


# Define the full parameter vector (initial parameters)
start_par <- c(
  0.9, 0,               # A row 1
  0, 0.9,               # A row 2
  0.1, 0.1, 0.1,              # B row 1
  0.1, 0.1, 0.1,              # B row 2
  0.1, 0, 0.1,          # Qlt (3 params)
  log(1)               # R_obs (log to ensure positive variance)
)


# Run estimation for the 2D model
result_2d <- estimate_dt(
  start_par = start_par,
  df = transformer_data,
  # Bounds: ensure diagonal of A stays within stable range (-1, 1)
  # and variances (Qlt, R_obs) stay positive
  lower = c(rep(-0.99, 4), rep(-Inf, 6), 1e-4, -Inf, 1e-4, -Inf),
  upper = c(rep(0.99, 4),  rep(Inf, 6),  Inf,  Inf,  Inf,  Inf)
)

# Extract the estimated model parameters
# The final parameter vector (par)
opt_par <- result_2d$par

# Reconstruct the matrices (the estimates)
A_est <- matrix(opt_par[1:4], nrow = 2, byrow = TRUE)
B_est <- matrix(opt_par[5:10], nrow = 2, byrow = TRUE)

# System and observation noise
Qlt_est <- matrix(0, 2, 2)
Qlt_est[1,1] <- opt_par[11]
Qlt_est[2,1] <- opt_par[12]
Qlt_est[2,2] <- opt_par[13]
Q_est <- Qlt_est %*% t(Qlt_est)
R_est <- exp(opt_par[14])


print(A_est)
print(B_est)
print(Qlt_est)

# Checking convergence to see if the initial values where okay
result_2d$convergence

# Checking the eigenvalues for A
eigen(A_est)$values
