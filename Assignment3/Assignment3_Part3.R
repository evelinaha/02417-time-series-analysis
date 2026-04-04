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
acf(values$Ph, 
    lag.max = 50, 
    main = "Autocorrelation of Ph", 
    na.action = na.pass)

# Cross-correlation function
ccf(values$Ph, values$Gv,
    lag.max = 50,
    main = "Cross-Correlation: Ph vs Gv", 
    na.action = na.pass)

# Cross-correlation function
ccf(values$Ph, values$Tdelta,
    lag.max = 50,
    main = "Cross-Correlation: Ph vs Tdelta", 
    na.action = na.pass)
