getwd()

library(ggplot2)
library(tidyr)
library(dplyr)

# 0. Read data
box_data <- read.csv("/Users/evelinahaite/Desktop/DTU_Masters/Time Series Analysis/Assignment3/assignment3_2026/box_data_60min.csv")

class(box_data$tdate)

head(box_data)

# assignment 3 part 3 

# Convert to datetime (correct format)
box_data$tdate <- as.POSIXct(box_data$tdate,
                             format = "%Y-%m-%d %H:%M:%S",
                             tz = "UTC")

# sanity check
box_data$tdate
class(box_data$tdate)

## Make the  variable a floating point (i.e.\ decimal number)
box_data$Ph <- as.numeric(box_data$Ph)
box_data$Tdelta <- as.numeric(box_data$Tdelta)
box_data$Gv <- as.numeric(box_data$Gv)

# part 3.1 plotting the data 

# PLOTTING THE 3 box_data IN A COMBINED PLOT 
par(mfrow = c(3,1))  # 3 rows, 1 column

plot(box_data$tdate, box_data$Ph, type="l", col="blue",
     main="Ph", xlab="Time", ylab="Ph")

plot(box_data$tdate, box_data$Tdelta, type="l", col="red",
     main="Tdelta", xlab="Time", ylab="Tdelta")

plot(box_data$tdate, box_data$Gv, type="l", col="green",
     main="Gv", xlab="Time", ylab="Gv")

#PART 3.2 TRAIN AND TESTING 
## Divide intro train and test set
testend <- as.POSIXct("2013-02-06 00:00:00", tz="UTC")
Data_train <- box_data[box_data$tdate <= testend, ]
Data_test <- box_data[box_data$tdate > testend, ]


cat("Length of Data_train:", nrow(Data_train), "\n")

# 3.3 Plotting scatter and correlations 

# 1. Autocorrelation of Ph
acf(Data_train$Ph, main = "Autocorrelation of Ph", lag.max = 50)

# 2. Cross-Correlation: Ph vs Gv
# Note: ccf(x, y) plots correlation between x[t+k] and y[t]
ccf(Data_train$Ph, Data_train$Gv, main = "Cross-Correlation: Ph vs Gv", lag.max = 50)

# 3. Cross-Correlation: Ph vs Tdelta
ccf(Data_train$Ph, Data_train$Tdelta, main = "Cross-Correlation: Ph vs Tdelta", lag.max = 50)


# 3.4 Impulse response
# Fit a regression model to see which lag values from the 2 independent variables are the most meaningful
model <- lm(Ph ~ Tdelta.l0 + Tdelta.l1 + Tdelta.l2 + Tdelta.l3 + Tdelta.l4 + Tdelta.l5 + 
              Tdelta.l6 + Tdelta.l7 + Tdelta.l8 + Tdelta.l9 + Tdelta.l10 +
              Gv.l0 + Gv.l1 + Gv.l2 + Gv.l3 + Gv.l4 + Gv.l5 + 
              Gv.l6 + Gv.l7 + Gv.l8 + Gv.l9 + Gv.l10 - 1, data = Data_train)

# Extract coefficients
results <- data.frame(term = names(coef(model)), estimate = as.numeric(coef(model)))

results

# Separate the variable name from the lag number
results <- results %>%
  separate(term, into = c("Variable", "Lag"), sep = "\\.l") %>%
  mutate(Lag = as.numeric(Lag))


ggplot(results, aes(x = Lag, y = estimate)) +
  geom_col(fill = "gray70", color = "black", width = 0.7) +
  facet_wrap(~Variable, scales = "free_y") +
  theme_bw() + 
  # Set breaks to be integers only
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + 
  labs(title = expression(paste("Estimated Impulse Response Weights (", psi[k], ")")),
       subtitle = "Response of Ph to Tdelta and Gv (Lags 0-10)",
       x = "Lag (k)",
       y = "Estimated Coefficient") +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        text = element_text(family = "serif"))


# 3.5 Multivariate linear model
# Fitting
multivar_model <- lm(Ph ~ Gv + Tdelta,
               data = Data_train)

summary(multivar_model)

# Print the coefficients
all_coefs <- coef(multivar_model)
all_coefs  

# Extract Residuals
Data_train$res <- residuals(multivar_model)

# 3. One-step prediction (fitted values)
Data_train$pred <- predict(multivar_model)

# 4. Plots
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Plot 1: Observed vs Predicted
plot(Data_train$Ph, type = "l", col = "gray70", lwd = 1,
     main = "One-step Prediction", ylab = "Ph (Heating Power)", xlab = "Time Index")
lines(Data_train$pred, col = "blue", lty = 2, lwd = 1.5)

legend("topleft", legend = c("Observed", "Predicted"), 
       col = c("gray70", "blue"), lty = 1:2, lwd = c(1, 1.5), 
       cex = 0.8) 

# Plot 2: Residuals over time
plot(Data_train$res, type = "l", col = "black",
     main = "Model Residuals", ylab = "Error (e_t)", xlab = "Time Index")
abline(h = 0, col = "red", lty = 3)



par(mfrow = c(2, 1), mar = c(4, 4, 4, 1))

acf(Data_train$res, main = "ACF of Residuals", lty = 2, lwd = 1.5, 
    xlab="Lag (hours)")
ccf(Data_train$res, Data_train$Tdelta, main = "CCF: Residuals vs Tdelta", lty = 2, lwd = 1.5, xlab="Lag (hours)", ylab = "CCF")


model_orders <- 1:10
rmse_results <- numeric(length(model_orders))
bic_results <- numeric(length(model_orders))
aic_results <- numeric(length(model_orders))

for(p in model_orders) {
  # Create a formula string
  ar_lags <- paste0("Ph.l", 1:p, collapse = " + ")
  form <- as.formula(paste("Ph ~ Gv + Tdelta +", ar_lags))
  
  # Fit on training data
  fit <- lm(form, data = Data_train)
  
  # Predict on test data
  preds <- predict(fit, newdata = Data_test)
  
  # Calculate RMSE
  residuals <- Data_test$Ph - preds
  rmse_results[p] <- sqrt(mean(residuals^2, na.rm = TRUE))
  # Store AIC and BIC results for comparison across lags
  aic_results[p] <- AIC(fit)
  bic_results[p] <- BIC(fit)
}

# 1. Identify the indices of the minimum values
best_rmse_idx <- which.min(rmse_results)
best_aic_idx  <- which.min(aic_results)
best_bic_idx  <- which.min(bic_results)

# 2. Set up the plotting grid
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

# Plot RMSE: Color all black, then the best one purple
rmse_cols <- rep("black", length(model_orders))
rmse_cols[best_rmse_idx] <- "red"

plot(model_orders, rmse_results, type = "b", pch = 19, col = rmse_cols,
     main = "RMSE vs Model Order", xlab = "Order (p)", ylab = "RMSE")
# Add a point on top to ensure the 'best' color is vivid
points(model_orders[best_rmse_idx], rmse_results[best_rmse_idx], col = "red", pch = 19, cex = 1.2)

# Plot AIC: Color best one red
aic_cols <- rep("black", length(model_orders))
aic_cols[best_aic_idx] <- "red"

plot(model_orders, aic_results, type = "b", pch = 19, col = aic_cols,
     main = "AIC vs Model Order", xlab = "Order (p)", ylab = "AIC")
points(model_orders[best_aic_idx], aic_results[best_aic_idx], col = "red", pch = 19, cex = 1.2)

# Plot BIC: Color best one red
bic_cols <- rep("black", length(model_orders))
bic_cols[best_bic_idx] <- "red"

plot(model_orders, bic_results, type = "b", pch = 19, col = bic_cols,
     main = "BIC vs Model Order", xlab = "Order (p)", ylab = "BIC")
points(model_orders[best_bic_idx], bic_results[best_bic_idx], col = "red", pch = 19, cex = 1.2)
# Lowest RMSE at Lag = 4, AIC = 6 and BIC = 5


# 3.9. Multi-Step Simulation 

# Final model based on RMSE results
fit_best <- lm(Ph ~ Gv + Tdelta + Ph.l1 + Ph.l2 + Ph.l3 + Ph.l4, data = box_data)

# Start simulation at index 5
n <- nrow(box_data)
sim_Ph <- numeric(n)
sim_Ph[1:4] <- box_data$Ph[1:4] 

for(t in 5:n) {
  # Note: Gv and Tdelta are always 'real' data
  # Ph lags are the 'simulated' values we generated in previous loop steps
  current_step <- data.frame(
    Gv = box_data$Gv[t],
    Tdelta = box_data$Tdelta[t],
    Ph.l1 = sim_Ph[t-1],
    Ph.l2 = sim_Ph[t-2],
    Ph.l3 = sim_Ph[t-3],
    Ph.l4 = sim_Ph[t-4]
  )
  
  sim_Ph[t] <- predict(fit_best, newdata = current_step)
  aic_results_best <- AIC(fit_best)
  bic_results_best <- BIC(fit_best)
}



# Plot the actual VS simulated

par(mfrow = c(1, 1))
plot(box_data$Ph,  type = "l", col = "black", alpha=0.7, lwd = 1,
     main = "Observed vs. Simulated in Multi-Step Simulation",
     xlab = "Time Index", ylab = "Heating Power (Ph)",
     ylim = range(c(box_data$Ph, sim_Ph), na.rm = TRUE)) 
lines(sim_Ph, col = "red", alpha=0.5, lwd = 1.5, lty = 1)


legend("topright", 
       legend = c("Observed (Actual)", "Multi-step Simulation"),
       col = c("black", "red", alpha=0.7), 
       lwd = c(1, 1), 
       bg = "white",
       cex = 0.8)



# Extract Residuals
box_data$res_multi <- residuals(sim_Ph)

# 3. One-step prediction (fitted values)
box_data$pred_multi <- predict(sim_Ph)

# 4. Plots
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Plot 1: Observed vs Predicted
plot(Data_train$Ph, type = "l", col = "gray70", lwd = 1,
     main = "One-step Prediction", ylab = "Ph (Heating Power)", xlab = "Time Index")
lines(Data_train$pred, col = "blue", lty = 2, lwd = 1.5)

legend("topleft", legend = c("Observed", "Predicted"), 
       col = c("gray70", "blue"), lty = 1:2, lwd = c(1, 1.5), 
       cex = 0.8) 

# Plot 2: Residuals over time
plot(Data_train$res, type = "l", col = "black",
     main = "Model Residuals", ylab = "Error (e_t)", xlab = "Time Index")
abline(h = 0, col = "red", lty = 3)



# One-step residuals (reusing the ones generated in section 3.5)
res_onestep <- Data_train$res

# Multi-step residuals
res_multistep <- box_data$Ph[5:n] - sim_Ph[5:n]


par(mfrow = c(2, 1), mar = c(4, 4, 4, 1))

# ACF of One-Step (Should look like white noise)
acf(res_onestep, main = "ACF of One-Step Residuals")

# ACF of Multi-Step (Will likely show high correlation)
acf(res_multistep, main = "ACF of Multi-Step Simulation Errors")
