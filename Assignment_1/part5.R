
D <- DST_BIL54
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

