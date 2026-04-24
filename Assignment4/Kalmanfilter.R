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
  P_filt  <- numeric(N)  # Filtered variances NOTE THIS IS ONLY VARIANCE HERE AS IT IS 1 DIMENSIONAL BUT BECOMES COVARIANCE IN 2D 
  innovation     <- numeric(N)  # Pre-fit residuals: y[t] - x_pred[t]
  innovation_var <- numeric(N)  # Innovation covariance: P_pred[t] + R
  
  for (t in seq_len(N)) {
    # the prediction step
    if (t == 1) {
      x_pred[t] <-  x_prior 
      P_pred[t] <-  P_prior
    } else {
      # calculations based on 10.63 in the book 
      x_pred[t] <- a %*% x_filt[t-1] + b %*% P_filt[t-1]  # the mean prediction using the previous filtered estimate
      P_pred[t] <- sigma1 %*% x_pred[t]                   # the variance prediction using the previous filtered estimate
    }
    
    # the update step
    innovation[t] <-  # the prediction error difference in measurement and predicted measurement 
    innovation_var[t] <- # the prediction error variance best guess :  St=Σt∣t1− +R
      
      
    K_t <-  P_filt[t-1] %*% t(sigma1) %*% solve(P_filt[t-1]) # the Kalman gain  10.75 in the book 
    x_filt[t] <- x_filt[t-1] + k_t %*% (y[t]-sigma1 %*% x_filt[t-1]) # the filtered estimate   10.73 
    P_filt[t] <- P_filt[t-1] - K_t %*% sigma1 %*% P_filt[t-1]   #   the filtered estimate variance 10.74 
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
