# =========================================================
# 1D Kalman Filter Log-Likelihood
# Same structure as the 2D version, but with:
#   - 1 state
#   - A: 1x1
#   - B: 1x3
#   - H: 1x1
#   - Q: scalar
# =========================================================

kf_logLik_dt <- function(par, df) {
  
  # -------------------------------------------------------
  # State-space matrices
  # -------------------------------------------------------
  
  # A is 1x1
  A <- matrix(par[1], nrow = 1, ncol = 1)
  
  # B is 1x3
  B <- matrix(par[2:4], nrow = 1, ncol = 3)
  
  # H is 1x1
  #H <- matrix(par[5], nrow = 1, ncol = 1)
  H <- matrix(1)
  
  obs_cols   <- c("Y")
  input_cols <- c("Ta", "S", "I")
  
  # -------------------------------------------------------
  # System noise covariance Q
  # Using exp() to ensure positivity
  # -------------------------------------------------------
  
  Q <- matrix(exp(par[6]), nrow = 1, ncol = 1)
  
  # -------------------------------------------------------
  # Observation noise covariance R
  # -------------------------------------------------------
  
  R_obs <- matrix(exp(par[7]), nrow = 1, ncol = 1)
  
  # -------------------------------------------------------
  # Initial state
  # -------------------------------------------------------
  
  X0 <- matrix(par[8], nrow = 1, ncol = 1)
  
  # -------------------------------------------------------
  # Data
  # -------------------------------------------------------
  
  Y  <- as.matrix(df[, obs_cols])
  U  <- as.matrix(df[, input_cols])
  Tn <- nrow(df)
  
  # -------------------------------------------------------
  # Initialization
  # -------------------------------------------------------
  
  n <- 1
  
  x_est <- X0
  P_est <- matrix(0.5, 1, 1)
  
  logLik <- 0
  
  
  # =======================================================
  # Kalman Filter Loop
  # =======================================================
  
  for (t in 1:Tn) {
    
    u_t <- matrix(U[t, ], 3, 1)
    
    # ---------------------------------------------------
    # Prediction step
    # ---------------------------------------------------
    
    if (t == 1) {
      
      x_pred <- X0
      P_pred <- matrix(10, 1, 1)
      
    } else {
      
      x_pred <- A %*% x_est + B %*% u_t
      P_pred <- A %*% P_est %*% t(A) + Q
      
    }
    
    # ---------------------------------------------------
    # Innovation
    # ---------------------------------------------------
    
    y_pred <- H %*% x_pred
    
    S_t <- H %*% P_pred %*% t(H) + R_obs
    
    innov <- matrix(Y[t], 1, 1) - y_pred
    
    # ---------------------------------------------------
    # Log-likelihood
    # ---------------------------------------------------
    
    logLik <- logLik - 0.5 * (
      log(det(S_t)) +
        t(innov) %*% solve(S_t) %*% innov
    )
    
    # ---------------------------------------------------
    # Update step
    # ---------------------------------------------------
    
    K_t <- P_pred %*% t(H) %*% solve(S_t)
    
    x_est <- x_pred + K_t %*% innov
    
    P_est <- (diag(1) - K_t %*% H) %*% P_pred
  }
  
  as.numeric(logLik)
}


# =========================================================
# Optimizer wrapper
# =========================================================

estimate_dt <- function(start_par, df, lower = NULL, upper = NULL) {
  
  negLL <- function(par) -kf_logLik_dt(par, df)
  
  optim(
    par     = start_par,
    fn      = negLL,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = 1000, trace = 1)
  )
}

# =========================================================
# Read data
# =========================================================

setwd("C:\\Users\\Frederik\\OneDrive - Danmarks Tekniske Universitet\\Uni\\Kandidat\\Time series Analysis")

transformer_data <- read.csv("transformer_data.csv")

# =========================================================
# Initial parameter vector
#
# par =
# [1]   A
# [2:4] B
# [5]   H
# [6]   log(Q)
# [7]   log(R)
# [8]   X0
# =========================================================

start_par <- c(
  0.9,                    # A
  1, 1, 1,                # B
  1,                      # H
  log(0.1),               # Q
  log(1),                 # R
  transformer_data$Y[1]   # X0
)

# =========================================================
# Run estimation
# =========================================================

result_1d <- estimate_dt(
  start_par = start_par,
  df = transformer_data,
  
  lower = c(
    -0.99,                # A
    rep(-Inf, 3),         # B
    -Inf,                 # H
    log(1e-4),            # Q
    log(1e-4),            # R
    -Inf                  # X0
  ),
  
  upper = c(
    0.99,                 # A
    rep(Inf, 3),          # B
    Inf,                  # H
    Inf,                  # Q
    Inf,                  # R
    Inf                   # X0
  )
)

# =========================================================
# Extract estimated parameters
# =========================================================

opt_par <- result_1d$par

A_est <- matrix(opt_par[1], 1, 1)

B_est <- matrix(opt_par[2:4], 1, 3)

H_est <- matrix(opt_par[5], 1, 1)

Q_est <- exp(opt_par[6])

R_est <- exp(opt_par[7])

X0_est <- matrix(opt_par[8], 1, 1)

# =========================================================
# Print results
# =========================================================

print(A_est)
print(B_est)
print(H_est)

print(Q_est)
print(R_est)

print(X0_est)

# =========================================================
# Convergence check
# 0 means successful convergence
# =========================================================

result_1d$convergence


# =========================================================
# Kalman Filter Function
# Obtain fitted values and residuals
# =========================================================

kalman_filter_run <- function(par, df) {
  
  # -------------------------------------------------------
  # Rebuild matrices
  # -------------------------------------------------------
  
  A <- matrix(par[1:4], nrow = 2, ncol = 2)
  B <- matrix(par[5:10], nrow = 2, ncol = 3)
  H <- matrix(par[11:12], nrow = 1, ncol = 2)
  
  Qlt <- matrix(0, 2, 2)
  Qlt[1,1] <- par[13]
  Qlt[2,1] <- par[14]
  Qlt[2,2] <- par[15]
  
  Q <- Qlt %*% t(Qlt)
  
  R_obs <- matrix(exp(par[16]))
  
  X0 <- matrix(par[17:18], nrow = 2, ncol = 1)
  
  # -------------------------------------------------------
  # Data
  # -------------------------------------------------------
  
  Y <- as.matrix(df$Y)
  U <- as.matrix(df[, c("Ta","S","I")])
  
  Tn <- nrow(df)
  
  # -------------------------------------------------------
  # Initialization
  # -------------------------------------------------------
  
  x_est <- X0
  P_est <- diag(1, 2)
  
  fitted_vals <- numeric(Tn)
  residuals   <- numeric(Tn)
  
  # =======================================================
  # Kalman loop
  # =======================================================
  
  for (t in 1:Tn) {
    
    u_t <- matrix(U[t, ], 3, 1)
    
    # Prediction
    if (t == 1) {
      
      x_pred <- X0
      P_pred <- diag(10, 2)
      
    } else {
      
      x_pred <- A %*% x_est + B %*% u_t
      P_pred <- A %*% P_est %*% t(A) + Q
      
    }
    
    # Predicted observation
    y_pred <- H %*% x_pred
    
    fitted_vals[t] <- y_pred[1,1]
    
    # Innovation
    innov <- Y[t] - y_pred
    
    residuals[t] <- innov[1,1]
    
    # Innovation covariance
    S_t <- H %*% P_pred %*% t(H) + R_obs
    
    # Kalman gain
    K_t <- P_pred %*% t(H) %*% solve(S_t)
    
    # Update
    x_est <- x_pred + K_t %*% innov
    
    P_est <- (diag(2) - K_t %*% H) %*% P_pred
  }
  
  list(
    fitted = fitted_vals,
    residuals = residuals
  )
}


# =========================================================
# Run filter with estimated parameters
# =========================================================

kf_results <- kalman_filter_run(opt_par, transformer_data)

residuals_kf <- kf_results$residuals
fitted_kf    <- kf_results$fitted

# =========================================================
# Information Criteria
# =========================================================

logLik_val <- kf_logLik_dt(opt_par, transformer_data)

k <- length(opt_par)
n <- nrow(transformer_data)

AIC_val <- -2 * logLik_val + 2 * k
BIC_val <- -2 * logLik_val + log(n) * k

cat("\n=============================\n")
cat("Model Selection Criteria\n")
cat("=============================\n")

cat("Log-Likelihood:", logLik_val, "\n")
cat("AIC:", AIC_val, "\n")
cat("BIC:", BIC_val, "\n")

# =========================================================
# Residual Diagnostics
# =========================================================

par(mfrow = c(2,2))

# Residual plot
plot(
  residuals_kf,
  type = "l",
  main = "Kalman Filter Residuals",
  ylab = "Residual",
  xlab = "Time"
)

abline(h = 0, col = "red")

# ACF
acf(
  residuals_kf,
  main = "ACF of Residuals"
)

# PACF
pacf(
  residuals_kf,
  main = "PACF of Residuals"
)

# QQ-plot
qqnorm(
  residuals_kf,
  main = "QQ-Plot of Residuals"
)

qqline(
  residuals_kf,
  col = "red"
)

# =========================================================
# Observed vs fitted
# =========================================================

par(mfrow = c(1,1))

plot(
  transformer_data$Y,
  type = "l",
  lwd = 2,
  main = "Observed vs Fitted",
  ylab = "Y",
  xlab = "Time"
)

lines(
  fitted_kf,
  col = "blue",
  lwd = 2
)

legend(
  "topright",
  legend = c("Observed", "Fitted"),
  col = c("black", "blue"),
  lty = 1,
  lwd = 2
)

