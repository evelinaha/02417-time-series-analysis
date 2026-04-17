# assignment 3 part 3 


## put in your own file dir. 
# Read data
setwd("C:\\Users\\Frederik\\OneDrive - Danmarks Tekniske Universitet\\Uni\\Kandidat\\Time series Analysis")
box_data <- read.csv("box_data_60min.csv")

# sanity check
str(box_data)

# Convert to datetime (correct format)
box_data$tdate <- as.POSIXct(box_data$tdate,
                           format = "%Y-%m-%d %H:%M:%S",
                           tz = "UTC")

# sanity check
box_data$tdate
class(box_data$tdate)

## Make the variable a floating point (i.e.\ decimal number), otherwise it's seen as a int. 
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

# debugging : 
#cat("Length of Data_train:", nrow(Data_train), "\n")

# 3.3 Plotting scatter and correlations 
plot(Data_test$Gv, Data_test$Ph,
     xlab = "Gv",
     ylab = "Ph",
     main = "Scatter Plot (Test Set)",
     pch = 16, col = "darkgreen")
abline(lm(Ph ~ Gv, data = values), col = "red", lwd = 2)

plot(Data_test$Tdelta, Data_test$Ph,
     xlab = "Ph",
     ylab = "Gv",
     main = "Scatter Plot (Test Set)",
     pch = 16, col = "darkgreen")
abline(lm(Ph ~ Tdelta, data = values), col = "red", lwd = 2)

# Basic autocorrelation plot of Ph
acf(Data_test$Ph, 
    lag.max = 50, 
    main = "Autocorrelation of Ph", 
    na.action = na.pass)

# Cross-correlation function
ccf(Data_test$Ph, Data_test$Gv,
    lag.max = 50,
    main = "Cross-Correlation: Ph vs Gv", 
    na.action = na.pass)

# Cross-correlation function
ccf(Data_test$Ph, Data_test$Tdelta,
    lag.max = 50,
    main = "Cross-Correlation: Ph vs Tdelta", 
    na.action = na.pass)



# 3.6 and 3.7  ARX model 
# Code provided by the solutions CANNOT MAKE IT WORK but confirms somewhat that what we are doing is correct 
lm(formula = ARX("Ph", c("Tdelta ", "Gv"), lags = 1),
   data = Data_train)
# Set the max order of the ARX modes that you wish to plot

# For easier understanding this was how the model looked like before being turned into a general while loop: 
arx1_lm <- lm(Ph ~ -Ph.l1 + Tdelta + Gv, data = Data_train)

arx7_lm <- lm(Ph ~ -(Ph.l1 + Ph.l2 + Ph.l3 + Ph.l4 + Ph.l5 + Ph.l6 + Ph.l7) + Tdelta + Tdelta.l1 + Tdelta.l2 + Tdelta.l3 + Tdelta.l4 + Tdelta.l5 + Tdelta.l6 + Gv + Gv.l1 + Gv.l2 + Gv.l3 + Gv.l4 + Gv.l5 + Gv.l6, data = Data_train)


# creating a one step prediction on the test data
One_step_predictions <- predict(arx1_lm, newdata = Data_test)

# plotting of the one step predictions: 
par(mfrow = c(1,1))
plot(Data_test$Ph, type = "l", col = "black", lwd = 2,
     ylab = "Ph [W]", xlab = "Time [Hours]", main = "One-step Predictions from first order ARX model on testdata")

lines(One_step_predictions, col = "red", lwd = 2)

legend("topleft",
       legend = c("Actual", "Predicted"),
       col = c("black", "red"),
       lty = 1,
       lwd = 2)


#debugging: 
#summary(arx1_lm)


# The ARX(p) model of order 2 
#arx2_lm <- lm(Ph ~ -(Ph.l1 + Ph.l2) + Tdelta + Tdelta.l1 + Gv + Gv.l1, data = Data_test) 
#summary(arx2_lm)




max_p <- 10
p <- 1

AIC_vals <- numeric(max_p)
BIC_vals <- numeric(max_p)
RMSE_results <- numeric(max_p)

# the while loop finds the ARX model for all ARX order up to p_max and then plots these. 
while (p <= max_p) {
  
  ar_terms <- paste0("Ph.l", 1:p, collapse = " + ")
  ar_terms <- paste0("-(", ar_terms, ")")
  
  tdelta_terms <- c("Tdelta", paste0("Tdelta.l", 1:(p-1)))
  tdelta_terms <- paste(tdelta_terms, collapse = " + ")
  
  gv_terms <- c("Gv", paste0("Gv.l", 1:(p-1)))
  gv_terms <- paste(gv_terms, collapse = " + ")
  
  formula_str <- paste("Ph ~", ar_terms, "+", tdelta_terms, "+", gv_terms)
  
  print(formula_str)
  
  model <- lm(as.formula(formula_str), data = Data_train)
  
  # Predict on test data
  preds <- predict(model, newdata = Data_test)
  
  # Calculate RMSE
  residuals <- Data_test$Ph - preds
  RMSE_results[p] <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  # calculate and store AIC and BIC scores 
  AIC_vals[p] <- AIC(model)
  BIC_vals[p] <- BIC(model)
  
  
  # loop counter
  p <- p + 1
  
  
}

par(mfrow = c(1,1))
# Plotting of the BIC and AIC values 
plot(1:max_p, AIC_vals, type = "b", pch = 19, col = "blue",
     ylim = range(c(AIC_vals, BIC_vals)),
     xlab = "Model Order (p)", 
     ylab = "Criterion Value",
     main = "AIC and BIC vs Model Order")
grid()
lines(1:max_p, BIC_vals, type = "b", pch = 17, lty = 2, col = "red")

legend("bottomleft",
       legend = c("AIC", "BIC"),
       col = c("blue", "red"),
       lty = c(1,2),
       pch = c(19,17))

# PLOTTING THE RESIDUALS of the First order ARX model  
par(mfrow = c(2,2))

# ACF: 
acf(residuals(arx1_lm), main = "ACF of Residuals")

# CCF between residuals and each predictor
ccf(residuals(arx1_lm), Data_test$Tdelta, main = "CCF: Residuals vs Tdelta")
ccf(residuals(arx1_lm), Data_test$Gv, main = "CCF: Residuals vs Gv")
ccf(residuals(arx1_lm), Data_test$Ph.l1, main = "CCF: Residals vs Ph.l1")


# PLOTTING THE RESIDUALS of the Fifth order ARX model (as is determined is the order with the lowest criterion score in the AIC and BIC plots)
par(mfrow = c(2,2))
# ACF: 
acf(residuals(arx1_lm), main = "ACF of Residuals")

# CCF between residuals and each predictor
ccf(residuals(arx1_lm), Data_test$Tdelta, main = "CCF: Residuals vs Tdelta")
ccf(residuals(arx1_lm), Data_test$Gv, main = "CCF: Residuals vs Gv")
ccf(residuals(arx1_lm), Data_test$Ph.l1, main = "CCF: Residals vs Ph.l1")

# 3.8 PLOTTING THE RMSE 

# 1. Identify the indices of the minimum values
best_rmse_idx <- which.min(RMSE_results)
best_aic_idx  <- which.min(AIC_vals)
best_bic_idx  <- which.min(BIC_vals)

# 2. Set up the plotting grid
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

# Plot RMSE: Color all black, then the best one purple
rmse_cols <- rep("black", length(max_p))
rmse_cols[best_rmse_idx] <- "red"


model_order <- 1:max_p

plot(model_order, RMSE_results, type = "b", pch = 19, col = rmse_cols,
     main = "RMSE vs Model Order", xlab = "Order (p)", ylab = "RMSE")
# Add a point on top to ensure the 'best' color is vivid
points(model_order[best_rmse_idx], RMSE_results[best_rmse_idx], col = "red", pch = 19, cex = 1.2)


# TESTING OF 3.9 

fit_best <- lm(Ph ~ -(Ph.l1 + Ph.l2 + Ph.l3 + Ph.l4) + Tdelta + Tdelta.l1 + Tdelta.l2 + Tdelta.l3 + Gv + Gv.l1 + Gv.l2 + Gv.l3, data = box_data)

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





