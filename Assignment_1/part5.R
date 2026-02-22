
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
plot(line_data_1$thetas[1, ], type='l', main='θ1,t lambda =  0.7')
plot(line_data_1$thetas[2, ], type='l', main='θ2,t, lambda =  0.7')
plot(line_data_2$thetas[1, ], type='l', main='θ1,t, lambda =  0.99')
plot(line_data_2$thetas[2, ], type='l', main='θ2,t, lambda =  0.99')

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
plot(residuals_1[5:length(residuals_1)], type='l', main='residuals with lambda =  0.7')
plot(residuals_2[5:length(residuals_2)], type='l', main='residuals with lambda =  0.99')

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
  return_value <- list("residuals" = residuals, "predictions" = predictions)
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
lambda = 0.6
R <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
RLS_model = RLS_with_forgetting(R, lambda)
final_theta = RLS_model$theta
