# R script for Part 4 of Assignment 1
# Imports
library(ggplot2)
### Read training data (stored in the global environment)
D <- read.csv("DST_BIL54.csv")
str(D)

D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

summary(D)

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]


# 4.1. Variance-covariance matrices: Global VS Local Model
# Setup for comparing matrices
N <- 72
lambda <- 0.9

# Weights for WLS (Local Model)
weights_local <- lambda^((N-1):0)

# Sigma for Local Model (WLS)
# Diagonal is 1/weight
Sigma_local <- diag(1 / weights_local)

# Sigma for Global Model (OLS)
Sigma_global <- diag(rep(1, N))

# The first 5 elements of the diagonal for both
cat("Global Sigma Diagonal (First 5):\n")
print(diag(Sigma_global)[1:5])

cat("\nLocal Sigma Diagonal (First 5):\n")
print(diag(Sigma_local)[1:5])


# 4.2. Lambda weights over time
plot(Dtrain$year, 1/diag(Sigma_local), 
     main='Lambda Weights VS Time (Training Data)',
     xlab='Time (Year)', ylab='Weight', col='hotpink', pch=16)
# Add a gridline for easier reading
grid()

# 4.3. Sum of all the lambda weights

weights <- 1/diag(Sigma_local)
sum_weights_wls <- sum(weights)


cat("The sum of lambda-weights for WLS is:", sum_weights_wls, "\n")

# Corresponding Sum of Weights in an OLS Model
sum_weights_ols <- sum(diag(Sigma_global))

cat("The sum of lambda-weights for OLS is:", sum_weights_ols, "\n")


# 4.4. Estimate and present parameters (theta)
X <- cbind(1, Dtrain$year)
SIGMA <- diag(weights)

# Get the predictor variable (total)
y <- Dtrain$total

# estimate parameters with WLS
theta_WLS <- solve(t(X) %*% SIGMA %*% X) %*% (t(X) %*% SIGMA %*% y)

cat("Local Intercept (theta1):", theta_WLS[1], "\n")
cat("Local Slope (theta2):", theta_WLS[2], "\n")


# FOR COMPARISON: Global OLS (Weight matrix is just Identity)
W_global <- diag(1, nrow = N) # Matrix of 1s on the diagonal
y <- Dtrain$total

# Solve Normal Equations for OLS
theta_global <- solve(t(X) %*% W_global %*% X) %*% (t(X) %*% W_global %*% y)

# Comparison
cat("GLOBAL OLS:\n")
print(theta_global)

cat("\nLOCAL WLS (lambda = 0.9):\n")
print(theta_WLS)


# 4.5. Forecast for the next 12 months
# Intercept and Slope from the WLS model
Tmemory <- sum(weights)
p <- 2
lambda <- 0.9
W <- diag(weights)

# WLS Estimation
XTWX_inv <- solve(t(X) %*% W %*% X)
theta_WLS <- XTWX_inv %*% t(X) %*% W %*% y

# Residual Variance
yhat_wls <- X %*% theta_WLS
e_wls <- y - yhat_wls
# Weighted RSS
RSS_wls <- sum(weights * e_wls^2) 
sigma2_wls <- as.numeric(RSS_wls / (Tmemory - p))

# Create test data placeholder for the next 12 months
last_time <- max(Dtrain$year)
future_months <- last_time + seq(1/12, 1, length.out=12)
Dtest_12 <- cbind(1, future_months)

# Point Forecasts and Variance
y_pred_wls <- Dtest_12 %*% theta_WLS

# Prediction Variance (includes sigma2_wls + parameter uncertainty)
V_pred <- sigma2_wls * (1 + diag(Dtest_12 %*% XTWX_inv %*% t(Dtest_12)))
se_pred <- sqrt(V_pred)

# Confidence Intervals (95%)
t_crit <- qt(0.975, df = Tmemory - p)
y_pred_lwr_wls <- y_pred_wls - t_crit * se_pred
y_pred_upr_wls <- y_pred_wls + t_crit * se_pred

# 6. Prepare data for vizualusation
df_forecast <- data.frame(
  year = future_months,
  y_pred = y_pred_wls,
  lwr = y_pred_lwr_wls,
  upr = y_pred_upr_wls,
  y_true = Dtest$total
)

print(head(df_forecast))


# OLS Comparison
# OLS Setup (Weights are all 1)
W_ols <- diag(1, nrow = N)

# OLS Estimation
XTX_inv <- solve(t(X) %*% X)
theta_OLS <- XTX_inv %*% t(X) %*% y

# Residual Variance
yhat_ols <- X %*% theta_OLS
e_ols <- y - yhat_ols
RSS_ols <- sum(e_ols^2) 
sigma2_ols <- as.numeric(RSS_ols / (N - p))

# Point Forecasts
y_pred_ols <- Dtest_12 %*% theta_OLS

# OLS Prediction Variance
V_pred_ols <- sigma2_ols * (1 + diag(Dtest_12 %*% XTX_inv %*% t(Dtest_12)))
se_pred_ols <- sqrt(V_pred_ols)

# Confidence Intervals (95%)
y_pred_lwr_ols <- y_pred_ols - t_crit * se_pred_ols
y_pred_upr_ols <- y_pred_ols + t_crit * se_pred_ols

# Prepare OLS forecast data for ggplot
df_forecast_ols <- data.frame(
  year = future_months,
  y_pred = y_pred_ols,
  lwr = y_pred_lwr_ols,
  upr = y_pred_upr_ols,
  y_true = Dtest$total
)

cat("\nLOCAL WLS (lambda = 0.9) Predictions:\n")
print(head(df_forecast))

cat("\nGLOBAL OLS Predictions:\n")
print(head(df_forecast_ols))


ggplot() +
  # True Data (Black)
  geom_point(data = df_forecast, aes(x = year, y = y_true, color = "True Data"), size = 1) +
  geom_line(data = df_forecast, aes(x = year, y = y_true, color = "True Data"), linetype = "dashed") +
  
  # OLS Forecast (Blue)
  geom_ribbon(data = df_forecast_ols, aes(x = year, ymin = lwr, ymax = upr, fill = "Global OLS"), alpha = 0.5) +
  geom_line(data = df_forecast_ols, aes(x = year, y = y_pred, color = "Global OLS"), size = 1.5) +
  
  # WLS Forecast (Pink)
  geom_ribbon(data = df_forecast, aes(x = year, ymin = lwr, ymax = upr, fill = "Local WLS"), alpha = 0.3) +
  geom_line(data = df_forecast, aes(x = year, y = y_pred, color = "Local WLS"), size = 1.5) +
  
  # True Data (Train) (Blue)
  geom_point(data = Dtrain, aes(x = year, y = total, color = "True Data"), size = 1) +
  geom_line(data = Dtrain, aes(x = year, y = total, color = "True Data"), linetype = "dashed") +
  
  # Define the Colors and Fills
  scale_color_manual(name = "Model Type", 
                     values = c("True Data" = "black", "Global OLS" = "cadetblue2", "Local WLS" = "hotpink")) +
  scale_fill_manual(name = "95% Interval", 
                    values = c("Global OLS" = "cadetblue2", "Local WLS" = "hotpink")) +
  
  # Labels
  labs(title = "Traffic Forecast Comparison: (OLS (Blue) VS WLS (Pink)",
       x = "Year", y = "Traffic (Millions)") +
  theme_minimal()

# Residuals
# Compute residuals for the test set
e_ols <- df_forecast_ols$y_true - df_forecast_ols$y_pred
e_wls <- df_forecast$y_true - df_forecast$y_pred

par(mfrow=c(1,2)) 
# Generate the Normal QQ-plot (OLS)
qqnorm(e_ols, main = "Q-Q Plot: OLS Residuals", col = "cadetblue2", pch = 16)
qqline(e_ols, col = "red", lwd = 2)

# Generate the WLS QQ-plot for comparison
qqnorm(e_wls, main = "Q-Q Plot: WLS Residuals", col = "hotpink", pch = 16)
qqline(e_wls, col = "blue", lwd = 2)

