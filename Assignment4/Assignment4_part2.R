

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







