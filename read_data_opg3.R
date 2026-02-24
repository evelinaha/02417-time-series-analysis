### Read training data



# initiating the dataframe 
?strftime
DST_BIL54$time <- as.POSIXct(paste0(DST_BIL54$time,"-01"), "%Y-%m-%d", tz="UTC")
DST_BIL54$time
class(DST_BIL54$time)

## Year to month for each of them
DST_BIL54$year <- 1900 + as.POSIXlt(DST_BIL54$time)$year + as.POSIXlt(DST_BIL54$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
DST_BIL54$total <- as.numeric(DST_BIL54$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- DST_BIL54[DST_BIL54$time < teststart, ]
Dtest <- DST_BIL54[DST_BIL54$time >= teststart, ]


# 1.1 Plot the training data
# Ensure we are using the fractional year (x-axis) and the total vehicles (y-axis)
par(mfrow = c(2, 2))

plot(Dtrain$year, Dtrain$total, 
     type = "o", 
     col = "hotpink",
     lwd = 2,
     main = " (Training Data)",
     cex.main = 0.9,    # Making the title (main) smaller
     xlab = "Time (Year)",
     ylab = "Total Vehicles (Millions)")
grid()


Dtrain$total[1:3]


# Part 3 
# the Y matrix is a 1x72 matrix with the values from the total column
Y <- Dtrain$total
# the X matrix is a 2X72 matrix with all 1's in the first column and the year timestamp in the 2nd column 
X <- model.matrix(~ year, data = Dtrain)

# We are here using the "manual" calculation method presented in the youtube
# video lectures 
F72 <- t(X)%*%X
h72 <- t(X)%*%Y

theta.hat72 <- solve(F72,h72)

## checking that the solution is correct with the R function lm()
model <- lm(total ~ year, data = Dtrain)


summary(model)

# initiating the intercept and slope parameter 
#Intercept 
theta1 <- theta.hat72[1,]
#Slope
theta2 <- theta.hat72[2,]

# calculating the standard error for each parameter using the YT lecture approach 
sigma2_hat <- as.numeric(t(Y-X %*% theta.hat72) %*% (Y-X %*% theta.hat72)/ 70)

# adding to the plot 
intercept <- theta.hat72[1]
slope     <- theta.hat72[2] 

abline(a = intercept, 
       b = slope, 
       col = "blue", 
       lwd = 2)

legend("topleft",
       legend = c("Observed data", "Fitted mean line"),
       col = c("hotpink", "blue"),
       lwd = 2,
       bty = "n")



# Predictions
# initiating the testing datapoints 
X_test <- model.matrix(~ year, data = Dtest)
# using the model parameters to make forecasted values from the test data
Y_predict <- as.numeric(X_test %*% theta.hat72)

#Confirming the predicted values with the R predict() function 
pred  <- predict(model, newdata = Dtest, interval = "prediction")

# PLOTTING THE FORECASTED VALUES 
# Combined y-limits
y_lim <- range(c(Dtest$total, Y_predict))

# Plot observed test data
plot(Dtest$year, Dtest$total,
     type = "o",
     col = "red",
     lwd = 2,
     pch = 16,
     ylim = y_lim,
     main = "(Testing Data: Observed vs Predicted)",
     cex.main = 0.9,
     xlab = "Time (Year)",
     ylab = "Total Vehicles (Millions)")
grid()

# Add predictions
lines(Dtest$year, Y_predict,
      col = "blue",
      lwd = 2,
      pch = 17,
      type = "o")

legend("topleft",
       legend = c("Observed (Test)", "Predicted"),
       col = c("red", "blue"),
       lwd = 2,
       pch = c(16, 17),
       bty = "n")


y_lim <- range(Dtest$total, pred[, "lwr"], pred[, "upr"])

plot(Dtest$year, Dtest$total,
     type = "p",
     col = "red",
     pch = 16,
     ylim = y_lim,
     main = "(Testing Data with Prediction Intervals)",
     cex.main = 0.9,
     xlab = "Time (Year)",
     ylab = "Total Vehicles (Millions)")
grid()

# Prediction interval bounds
lines(Dtest$year, pred[, "lwr"],
      col = "darkgreen",
      lwd = 2,
      lty = 2)

lines(Dtest$year, pred[, "upr"],
      col = "darkgreen",
      lwd = 2,
      lty = 2)
# Predicted mean
lines(Dtest$year, pred[, "fit"],
      col = "blue",
      lwd = 2)

legend("topleft",
       legend = c("Observed (Test)",
                  "Predicted Mean",
                  "95% Prediction Interval"),
       col = c("red", "blue", "darkgreen"),
       pch = c(16, NA, NA),
       lwd = c(NA, 2, 2),
       lty = c(NA, 1, 2),
       bty = "n")
# Polygon for shaded interval
polygon(c(Dtest$year, rev(Dtest$year)),
        c(pred[, "lwr"], rev(pred[, "upr"])),
        col = rgb(0, 0.6, 0, 0.2),
        border = NA)

# Redraw mean line on top
lines(Dtest$year, pred[, "fit"], col = "blue", lwd = 2)

residuals(model)
predict(model)
summary(model) 

#plot(model, which=1, col=c("blue")) # Residuals vs Fitted Plot


plot(model, which=2, col=c("red"))  # Q-Q Plot

