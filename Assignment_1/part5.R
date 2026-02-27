
D <- DST_BIL54
D <- read.csv("DST_BIL54.csv")
str(D)

# See the help
?strftime
D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")
D$time
class(D$time)

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]


# Part 5

# Initialize R and theta
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
theta <- matrix(c(0, 0), nrow = 2, ncol = 1)

# 5.2 Recursive Least Squares (RLS) implementation for the first 3 steps


# Create a matrix to store the theta estimates for the 3 time steps
theta_estimates <- matrix(NA, nrow = 2, ncol = 3)

# Run the loop
for (t in 1:3) {
  # Extract the current input x_t and output y_t
  x_t <- matrix(c(1, Dtrain$year[t]), nrow = 2, ncol = 1)
  y_t <- Dtrain$total[t]
  
  # Update R_t
  R <- R + x_t %*% t(x_t)
  
  # Update theta_t
  # solve(R) calculates the inverse of matrix R
  prediction_error <- y_t - (t(x_t) %*% theta)[1,1]   # Adding [1,1] at the end to extract the number of the matrix object (R stores it as a 1x1 matrix)
  theta <- theta + solve(R) %*% x_t * prediction_error
  
  # Store the estimate
  theta_estimates[, t] <- theta
}

# Print the estimates up to t=3
print("Theta estimates for t=1, 2, 3:")
print(theta_estimates)


# 5.3 Calculate RLS estimates for the entire training set (t = N)
N <- nrow(Dtrain)

# Re-initialize with the original values
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
theta <- matrix(c(0, 0), nrow = 2, ncol = 1)

# Run the loop for all N time steps
for (t in 1:N) {
  x_t <- matrix(c(1, Dtrain$year[t]), nrow = 2, ncol = 1)
  y_t <- Dtrain$total[t]
  
  R <- R + x_t %*% t(x_t)
  prediction_error <- y_t - (t(x_t) %*% theta)[1,1]
  theta <- theta + solve(R) %*% x_t * prediction_error
}

# Print the final RLS estimate
print("Final RLS estimate at t=N:")
print(theta)

#------------------------------------------------------------------------#

# 5.4
lambda_1 <- 0.7
lambda_2 <- 0.99

RLS_with_forgetting <- function (R, lambda) {
  theta <- matrix(c(0, 0), nrow = 2, ncol = 1)
  thetas <- matrix(0, nrow = 2, ncol = N)
  for (t in 1:N) {
    x_t <- matrix(c(1, Dtrain$year[t]), nrow = 2, ncol = 1)
    y_t <- Dtrain$total[t]
    
    R <- (lambda*R) + x_t %*% t(x_t)
    prediction_error <- y_t - (t(x_t) %*% theta)[1,1]
    theta <- theta + solve(R) %*% x_t * prediction_error
    thetas[, t] = theta
  }
  return_values = list("theta" = theta, "thetas" = thetas)
  return(return_values)
}

# Initialize R
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
line_data_1 <- RLS_with_forgetting(R, lambda_1)

# Initialize R
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
line_data_2 <- RLS_with_forgetting(R, lambda_2)

par(mfrow=c(2, 2))
plot(Dtrain$year, line_data_1$thetas[1, ], type='l', main='θ1,t lambda =  0.7', xlab='Year', ylab='Theta value')
plot(Dtrain$year, line_data_1$thetas[2, ], type='l', main='θ2,t, lambda =  0.7', xlab='Year', ylab='Theta value')
plot(Dtrain$year, line_data_2$thetas[1, ], type='l', main='θ1,t, lambda =  0.99', xlab='Year', ylab='Theta value')
plot(Dtrain$year, line_data_2$thetas[2, ], type='l', main='θ2,t, lambda =  0.99', xlab='Year', ylab='Theta value')

#------------------------------------------------------------------------#

# 5.5
one_step_predictions <- function (R, lambda) {
  theta <- matrix(c(0, 0), nrow = 2, ncol = 1)
  thetas <- matrix(0, nrow = 2, ncol = N)
  predictions <- matrix(NA, nrow=N)
  residuals <- matrix(NA, nrow=N)
  for (t in 1:N) {
    x_t <- matrix(c(1, Dtrain$year[t]), nrow = 2, ncol = 1)
    x_t_1 <- matrix(c(1, Dtrain$year[t+1]), nrow = 2, ncol = 1)
    y_t <- Dtrain$total[t]
    y_t_1 <- Dtrain$total[t+1]
    
    R <- (lambda*R) + x_t %*% t(x_t)
    prediction_error <- y_t - (t(x_t) %*% theta)[1,1]
    theta <- theta + solve(R) %*% x_t * prediction_error
    thetas[, t] = theta
    
    one_step_prediction <- (t(x_t_1)%*%theta)[1, 1]
    one_step_prediction_error <- y_t_1 - one_step_prediction
    residuals[t+1] = one_step_prediction_error
    predictions[t+1] = one_step_prediction
  }
  return_value <- list("residuals" = residuals, "predictions" = predictions)
  return(return_value)
}

# run function for both lambdas
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
prediction_data_1 = one_step_predictions(R, lambda_1)
predictions_1 = prediction_data_1$predictions
residuals_1 = prediction_data_1$residuals

R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
prediction_data_2 = one_step_predictions(R, lambda_2)
predictions_2 = prediction_data_2$predictions
residuals_2 = prediction_data_2$residuals

par(mfrow=c(2, 1))
plot(residuals_1[5:length(residuals_1)], type='l', main='residuals with lambda =  0.7', xlab='Observation', ylab='Residual')
plot(residuals_2[5:length(residuals_2)], type='l', main='residuals with lambda =  0.99', xlab='Observation', ylab='Residual')

#------------------------------------------------------------------------#

# 5.6
k_step_predictions <- function (R, lambda, k) {
  theta <- matrix(c(0, 0), nrow = 2, ncol = 1)
  thetas <- matrix(0, nrow = 2, ncol = N)
  predictions <- matrix(NA, nrow=N)
  residuals <- matrix(NA, nrow=N)
  for (t in 1:(N-k)) {
    x_t <- matrix(c(1, Dtrain$year[t]), nrow = 2, ncol = 1)
    y_t <- Dtrain$total[t]
    
    x_tk <- matrix(c(1, Dtrain$year[t+k]), nrow = 2, ncol = 1)
    y_tk <- Dtrain$total[t+k]
    
    R <- (lambda*R) + x_t %*% t(x_t)
    prediction_error <- y_t - (t(x_t) %*% theta)[1,1]
    theta <- theta + solve(R) %*% x_t * prediction_error
    thetas[, t] = theta
    
    k_step_prediction <- (t(x_tk)%*%theta)[1, 1]
    k_step_prediction_error <- y_tk - k_step_prediction
    residuals[t+k] = k_step_prediction_error
    predictions[t+k] = k_step_prediction
  }
  
  return_value <- list("residuals" = residuals, "predictions" = predictions, "theta" = theta)
  return(return_value)
}

# run function for all k and all lambda
NUMBER_OF_K <- 12
lambdas <- seq(0.5, 1.0, by = 0.01)
rmse_matrix <- matrix(NA, nrow = length(lambdas), ncol = NUMBER_OF_K)
for (k in 1:NUMBER_OF_K) {
  for (i in seq_along(lambdas)) {
    lambda <- lambdas[i]
    R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
    k_step_prediction_data <- k_step_predictions(R, lambda, k)
    k_residuals <- k_step_prediction_data$residuals
    RMSEk = sqrt(mean(k_residuals^2, na.rm = TRUE))
    rmse_matrix[i, k] <- RMSEk
  }
}

# Google Gemini used to make the plot
# --- Plotting ---

# 1. Define colors for the 12 different lines
colors <- rainbow(NUMBER_OF_K)

# 2. Initialize the plot using the first k (k=1)
# we set ylim to cover the full range of all data
plot(lambdas, rmse_matrix[, 1], 
     type = "l", 
     col = colors[1], 
     ylim = range(rmse_matrix, na.rm = TRUE),
     xlab = "Lambda", 
     ylab = "RMSEk", 
     main = "RMSE vs Lambda for various k-steps")

# 3. Add lines for k = 1 through 12
for (k in 1:NUMBER_OF_K) {
  lines(lambdas, rmse_matrix[, k], col = colors[k])
}

# 4. Add a legend to identify the lines
legend("topleft", 
       legend = paste("k =", 1:NUMBER_OF_K), 
       col = colors, 
       lty = 1, 
       cex = 0.7, 
       ncol = 2) # ncol=2 makes the legend more compact

#------------------------------------------------------------------------#

# 5.7
lambda = 0.7
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
RLS = RLS_with_forgetting(R, lambda)
final_theta = RLS$theta
predicted_values <- numeric(12)


for (t in 1:12) {
  x_t <- matrix(c(1, Dtest$year[t]), nrow = 2, ncol = 1)
  predicted_values[t] <- (t(x_t) %*% final_theta)[1,1]
}

plot(Dtest$year, predicted_values, 
     type = "l", 
     col = "blue", 
     lwd = 2,
     main = "Predicted vs Actual Values (Test Set)",
     xlab = "Time",
     ylab = "Value")

lines(Dtest$year, Dtest$total,
      col = "red",
      lwd = 2)

legend("topleft",
       legend = c("Predicted", "Actual"),
       col = c("blue", "red"),
       lwd = 2)

# ============================================================
# COMPARISON: OLS vs WLS vs RLS on the Test Set (Generated with Claude)
# ============================================================

# --- RLS final theta (from file 2, lambda = 0.7) ---
lambda_rls <- 0.9
R_init <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
rls_result <- RLS_with_forgetting(R_init, lambda_rls)
theta_RLS <- rls_result$theta

# --- Generate predictions for each model on Dtest ---
X_test <- cbind(1, Dtest$year)

pred_ols <- X_test %*% theta_OLS
pred_wls <- X_test %*% theta_WLS
pred_rls <- X_test %*% theta_RLS

# --- RMSE Comparison ---
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

cat("Test Set RMSE Comparison:\n")
cat("  OLS RMSE:", round(rmse(Dtest$total, pred_ols), 4), "\n")
cat("  WLS RMSE:", round(rmse(Dtest$total, pred_wls), 4), "\n")
cat("  RLS RMSE:", round(rmse(Dtest$total, pred_rls), 4), "\n")

# --- Combined Plot ---
y_range <- range(c(Dtest$total, pred_ols, pred_wls, pred_rls))

plot(Dtest$year, Dtest$total,
     type = "b", pch = 16, col = "black", lwd = 2,
     ylim = y_range,
     main = "Test Set Forecast: OLS vs WLS vs RLS",
     xlab = "Year", ylab = "Traffic (Millions)")

lines(Dtest$year, pred_ols, col = "cadetblue2", lwd = 2, lty = 2)
lines(Dtest$year, pred_wls, col = "hotpink",    lwd = 2, lty = 2)
lines(Dtest$year, pred_rls, col = "darkorange",  lwd = 2, lty = 2)

legend("topleft",
       legend = c("Actual", 
                  paste0("Global OLS"), 
                  paste0("Local WLS (λ=0.9)"),
                  paste0("RLS (λ=", lambda_rls, ")")),
       col    = c("black", "cadetblue2", "hotpink", "darkorange"),
       lwd    = 2, lty = c(1, 2, 2, 2), pch = c(16, NA, NA, NA))


