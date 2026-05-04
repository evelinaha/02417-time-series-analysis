

# assignment 4 part 2 
#install.packages("ggplot2")
#install.packages("tidyr")

## put in your own file dir. 
# Read data
setwd("C:\\Users\\Frederik\\OneDrive - Danmarks Tekniske Universitet\\Uni\\Kandidat\\Time series Analysis")
transformer_data <- read.csv("transformer_data.csv")

# sanity check
str(transformer_data)

# sanity check
transformer_data$time
class(transformer_data$time)

## Make the variable a floating point (i.e.\ decimal number), otherwise it's seen as a int. 
transformer_data$Y <- as.numeric(transformer_data$Y)
transformer_data$Ta <- as.numeric(transformer_data$Ta)
transformer_data$S <- as.numeric(transformer_data$S)
transformer_data$I <- as.numeric(transformer_data$I)


#################### PLOTTING THE DATA 4.1 #####################################

# Combine variables into a matrix
data_matrix <- cbind(transformer_data$Y,
                     transformer_data$Ta,
                     transformer_data$S,
                     transformer_data$I)

# Plot
matplot(transformer_data$time, data_matrix,
        type = "l", lty = 1, lwd = 2,
        col = c("blue", "red", "green", "purple"),
        xlab = "Time", ylab = "Values")

# Add legend
legend("topright",
       legend = c("Y", "Ta", "S", "I"),
       col = c("blue", "red", "green", "purple"),
       lty = 1, lwd = 2)
# Plot Y, Ta, I on left axis
matplot(transformer_data$time,
        cbind(transformer_data$Y,
              transformer_data$Ta,
              transformer_data$I),
        type = "l", lty = 1, lwd = 2,
        col = c("blue", "red", "purple"),
        xlab = "Time", ylab = "Y / Ta / I")

# Add S on right axis
par(new = TRUE)

plot(transformer_data$time,
     transformer_data$S,
     type = "l", lwd = 2,
     col = "green",
     axes = FALSE, xlab = "", ylab = "")

# Right axis
axis(side = 4)
mtext("S values", side = 4, line = 3)

# Left legend (for Y, Ta, I)
legend("topleft",
       legend = c("Y", "Ta", "I"),
       col = c("blue", "red", "purple"),
       lty = 1, lwd = 2,
       bty = "n")

# Right legend (for S)
legend("topright",
       legend = "S",
       col = "green",
       lty = 1, lwd = 2,
       bty = "n")


#######################  

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
  
  
  # debugging 
  print( dim(A))
  print(dim(x_est))
  print(dim(B))
  print(dim(U[t, ]))
  
  
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
      # the prediction step
      if (t == 1) {
        x_pred <- X0
        P_pred <- P_est
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



## --- Select relevant columns (optional but clean) ---
df <- transformer_data[, c("Y", "Ta", "S", "I")]

## --- Define starting parameters ---
# Example for 2D state:
# A (2x2) = 4 params
# B (2x3) = 6 params
# Qlt (lower triangular 2x2) = 3 params
# R_obs (1x1) = 1 param
start_par <- c(
  0.5, 0,   # A row 1
  0, 0.5,   # A row 2
  0.1, 0.1, 0.1,   # B row 1
  0.1, 0.1, 0.1,   # B row 2
  0.1, 0, 0.1,     # Qlt (lower triangular)
  0.1              # R_obs
)

## --- Test log-likelihood ---
ll_val <- kf_logLik_dt(start_par, df)
print(ll_val)

## --- Run estimation ---
result <- estimate_dt(
  start_par = start_par,
  df = df
)

## --- Output results ---
print(result$par)
print(result$value)



