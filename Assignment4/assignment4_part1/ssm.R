#### 1.1
set.seed(123)

a <- 0.9
b <- 1
sigma1 <- 1
X0 <- 5
n <- 100
num_of_realisations <- 5

# Matrix to store trajectories
# rows = time 0..100, columns = realization 1..5
X <- matrix(0, nrow = n + 1, ncol = num_of_realisations)

# Initial value
X[1, ] <- X0

# Simulate 5 independent realizations
for (j in 1:num_of_realisations) {
  for (t in 2:(n + 1)) {
    e_t <- rnorm(1, mean = 0, sd = sigma1)
    X[t, j] <- a * X[t - 1, j] + b + e_t
  }
}

#----------------------------------------------#
### PLOT
time <- 0:n
plot(time, X[,1],
     type = "l",
     col = 1,
     lwd = 2,
     ylim = range(X),
     xlab = "Time",
     ylab = "X_t",
     main = "5 Independent Realizations of State Process")

# Add remaining trajectories
for (j in 2:num_of_realisations) {
  lines(time, X[,j], col = j, lwd = 2)
}

# Legend
legend("topleft",
       legend = paste("Path", 1:num_of_realisations),
       col = 1:num_of_realisations,
       lwd = 2)

#----------------------------------------------#

### 1.2
sigma2 <- 1

X <- numeric(n + 1)
Y <- numeric(n + 1)

X[1] <- X0

Y[1] <- X[1] + rnorm(1, mean = 0, sd = sigma2)

for (t in 2:(n + 1)) {
  e1_t <- rnorm(1, mean = 0, sd = sigma1)
  X[t] <- a * X[t - 1] + b + e1_t
  
  e2_t <- rnorm(1, mean = 0, sd = sigma2)
  Y[t] <- X[t] + e2_t
}

#----------------------------------------------#
### PLOT
time <- 0:n

plot(time, X,
     type = "l",
     col = "blue",
     pch = 16,
     lwd = 2,
     ylim = range(c(X, Y)),
     xlab = "Time",
     ylab = "Value",
     main = "Latent State X_t and Observations Y_t")

lines(time, Y,
      type = "l",
      col = "red",
      pch = 16,
      lwd = 2)

legend("topleft",
       legend = c("Latent state X_t", "Observations Y_t"),
       col = c("blue", "red"),
       lty = 1,
       lwd = 2)

#----------------------------------------------#

### 1.3
# Kalman filter
myKalmanFilter <- function(
    y,             # Vector of observations y_t
    theta,       # Model parameters for X_{t+1} = a - b*X_t + c*e_t
    R,             # Measurement noise variance
    x_prior = 0,   # Initial prior mean for X_0
    P_prior = 10   # Initial prior variance for X_0
) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  N <- length(y)
  x_pred  <- numeric(N)  # Predicted means
  P_pred  <- numeric(N)  # Predicted variances
  x_filt  <- numeric(N)  # Filtered means
  P_filt  <- numeric(N)  # Filtered variances
  innovation     <- numeric(N)  # Pre-fit residuals: y[t] - x_pred[t]
  innovation_var <- numeric(N)  # Innovation covariance: P_pred[t] + R
  
  for (t in seq_len(N)) {
    # the prediction step
    if (t == 1) {
      x_pred[t] <- x_prior
      P_pred[t] <- P_prior
    } else {
      x_pred[t] <- a*x_filt[t-1] + b
      P_pred[t] <- a*P_filt[t-1]*t(a) + sigma1 # Maybe it shouldn't be power of 2
    }
    
    # the update step
    innovation[t] <- y[t] - x_pred[t]
    innovation_var[t] <- P_pred[t] + R
    
    K_t <- P_pred[t] * solve(innovation_var[t])
    
    x_filt[t] <- x_pred[t] + K_t * (y[t]-x_pred[t])
    P_filt[t] <- P_pred[t] - K_t*innovation_var[t]*t(K_t)
  }
  
  return(list(
    x_pred = x_pred,
    P_pred = P_pred,
    x_filt = x_filt,
    P_filt = P_filt,
    innovation = innovation,
    innovation_var = innovation_var
  ))
}

# Run Kalman filter
theta <- c(a, b, sigma1)
R <- sigma2^2

kf <- myKalmanFilter(
  y = Y,
  theta = theta,
  R = R,
  x_prior = X0,
  P_prior = 10
)

# 95% confidence interval around predicted state
lower <- kf$x_pred - 1.96 * sqrt(kf$P_pred)
upper <- kf$x_pred + 1.96 * sqrt(kf$P_pred)

#----------------------------------------------#
### PLOT
# Plot everything
plot(time, X,
     type = "b",
     col = "blue",
     pch = 16,
     cex = 0.6,
     lwd = 2,
     ylim = range(c(X, Y, lower, upper)),
     xlab = "Time",
     ylab = "Value",
     main = "Kalman Filter: State, Observations, Prediction, and 95% CI")

# Observation Y_t
lines(time, Y,
      type = "b",
      col = "red",
      pch = 16,
      cex = 0.6,
      lwd = 1)

# Predicted state x_hat_{t|t-1}
lines(time, kf$x_pred,
      type = "l",
      col = "darkgreen",
      lwd = 2,
      lty = 2)

# 95% confidence interval
lines(time, lower,
      col = "darkgreen",
      lwd = 1,
      lty = 3)

lines(time, upper,
      col = "darkgreen",
      lwd = 1,
      lty = 3)

legend("topleft",
       legend = c("Latent state X_t",
                  "Observation Y_t",
                  "Predicted state",
                  "95% CI"),
       col = c("blue", "red", "darkgreen", "darkgreen"),
       lty = c(1, 1, 2, 3),
       pch = c(16, 16, NA, NA),
       lwd = c(2, 1, 2, 1))

#----------------------------------------------#

### 1.4
myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  
  kf_result <- myKalmanFiter(y, theta, R, x_prior, P_prior)
  err <- kf_result$innovation
  S <- kf_result$innovatizon_var
  
  # Compute log-likelihood contributions from each time step
  logL <- sum(dnorm(err,
                    mean = 0,
                    sd = sqrt(S),
                    log = TRUE))
  logL <- 
    return(-logL)  # Return negative log-likelihood for minimization
}

# Known observation noise variance
R <- 1

# Initial guess: (a, b, sigma1)
theta_start <- c(a, b, sigma1)

# Estimate parameters by minimizing negative log-likelihood
fit <- optim(
  par = theta_start,
  fn = myLogLikFun,
  y = Y,
  R = R,
  x_prior = X0,
  P_prior = 10,
  method = "L-BFGS-B",
  lower = c(-Inf, -Inf, 0.001),   # sigma1 must be positive
  upper = c(Inf, Inf, Inf)
)

# Print results
cat("Estimated a      =", fit$par[1], "\n")
cat("Estimated b      =", fit$par[2], "\n")
cat("Estimated sigma1 =", fit$par[3], "\n")

############################
# Simulate many realizations for n=100
############################
new_a <- 1
new_b <- 0.9
new_sigma1 <- 1

X_mat <- matrix(NA, nrow = n + 1, ncol = n)
Y_mat <- matrix(NA, nrow = n + 1, ncol = n)

X_mat[1, ] <- X0

for (j in 1:n) {
  
  Y_mat[1, j] <- X_mat[1, j] + rnorm(1, 0, R)
  
  for (t in 2:(n + 1)) {
    
    X_mat[t, j] <- new_a * X_mat[t - 1, j] + new_b + rnorm(1, 0, new_sigma1)
    
    Y_mat[t, j] <- X_mat[t, j] + rnorm(1, 0, sigma2)
  }
}

############################
# Estimate parameters
############################
estimates <- matrix(NA, nrow = n, ncol = 3)
colnames(estimates) <- c("a_hat", "b_hat", "sigma1_hat")

theta_start <- c(0.8, 0.8, 1)

for (j in 1:n) {
  
  fit <- optim(
    par = theta_start,
    fn = myLogLikFun,
    y = Y_mat[, j],
    R = R,
    x_prior = X0,
    P_prior = 10,
    method = "L-BFGS-B",
    lower = c(-Inf, -Inf, 0.001),
    upper = c(Inf, Inf, Inf)
  )
  
  estimates[j, ] <- fit$par
}

############################
# Numerical summary
############################
cat("Means of estimates:\n")
print(colMeans(estimates))

cat("\nStandard deviations:\n")
print(apply(estimates, 2, sd))

cat("\nFive-number summaries:\n")
print(summary(estimates))

############################
# Boxplots
############################
par(mfrow = c(1, 3))

boxplot(estimates[, "a_hat"],
        main = "Estimate of a",
        ylab = expression(hat(a)),
        col = "lightblue")
abline(h = new_a, col = "red", lty = 2, lwd = 2)

boxplot(estimates[, "b_hat"],
        main = "Estimate of b",
        ylab = expression(hat(b)),
        col = "lightgreen")
abline(h = new_b, col = "red", lty = 2, lwd = 2)

boxplot(estimates[, "sigma1_hat"],
        main = expression("Estimate of " * sigma[1]),
        ylab = expression(hat(sigma)[1]),
        col = "lightpink")
abline(h = new_sigma1, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))

### 1.5

### 1.5 bonus
