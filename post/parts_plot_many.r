# File:   parts_plot_many.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the mean-square-displacement of the healthy (blue) and cancerous (red) cells by 
#   averaging multiple trials. A log-log slope is displayed with uncertainty.

# Read in data
base <- "~/Documents/CancerModeling/outputs/testB"
data_location <- "~/Documents/CancerModeling/outputs/testB1.csv"
sim_data <- read.csv(data_location, comment.char="#", sep="\t", header=FALSE)
time <- unique(sim_data$V1)

# Setup time and MSD vectors
healthy_msd <- numeric(length(time))
cancer_msd <- numeric(length(time))

# start full loop
num_trials <- 30
slope_H <- numeric(num_trials)
slope_C <- numeric(num_trials)
for (f in seq.int(1,num_trials)) {
  name <- paste(f,".csv", sep="")
  data_location <- paste(base,name,sep="")
  sim_data <- read.csv(data_location, comment.char="#", sep="\t", header=FALSE)
  time <- unique(sim_data$V1)
  
  # Get initial locations
  data <- subset(sim_data, V1 == time[1])
  healthy0 <- subset(data, V3 == "H")[,4:6]
  cancer0 <- subset(data, V3 == "C")[,4:6]

  i <- 1
  for (t in time) {
    data <- subset(sim_data, V1 == t)
    healthy <- subset(data, V3 == "H")
    cancer <- subset(data, V3 == "C")
    healthy_msd[i] <- mean(mapply(function(x, x0, y, y0, z, z0){(x-x0)^2+(y-y0)^2+(z-z0)^2},
                               healthy$V4, healthy0[,1], healthy$V5, healthy0[,2], 
                               healthy$V6, healthy0[,3] ))
    cancer_msd[i] <- mean(mapply(function(x, x0, y, y0, z, z0){(x-x0)^2+(y-y0)^2+(z-z0)^2},
                                cancer$V4, cancer0[,1], cancer$V5, cancer0[,2], 
                                cancer$V6, cancer0[,3] ))
    i <- i + 1
  }
  
  
  
  # Linear Regression
  #Start log fits at t=dt
  l_time = log(time[2:length(time)])
  l_msd_H = log(healthy_msd[2:length(time)])
  l_msd_C = log(cancer_msd[2:length(time)])
  
  # Linear Regression
  hfit <- lm(l_msd_H ~ l_time)
  cfit <- lm(l_msd_C ~ l_time)
  
  slope_H[f] <- hfit$coefficients[2]
  slope_C[f] <- cfit$coefficients[2]
}

# Plotting data with Linear Regression fitting
plot(l_time, l_msd_H, pch = 16, cex = 1.3, col = "blue", main = "HEALTHY-CELL MSD",
     xlab = "Time (tau)", ylab = "MSD (sigma^2)")
abline(hfit)
plot(l_time, l_msd_C, pch = 16, cex = 1.3, col = "red", main = "CANCER-CELL MSD")
abline(cfit)


print(sprintf("Healthy Slope: %.3f", hfit$coefficients[2]))
print(sprintf("Cancer Slope: %.3f", cfit$coefficients[2]))


print(sprintf("Healthy Slope: %.3f +/- %.3f", mean(slope_H), sd(slope_H)))
print(sprintf("Cancer Slope: %.3f +/- %.3f", mean(slope_C), sd(slope_C)))