# File:   calc_dfc.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the distance-from-center of the healthy (blue) and cancerous (red) cells.

# Read in data
data_location <- "~/Documents/CancerModeling/outputs/testA.csv"
sim_data <- read.csv(data_location, comment.char="#", sep="\t", header=FALSE)

# Setup time and MSD vectors
time <- unique(sim_data$V1)
healthy_dfc <- numeric(length(time))
cancer_dfc <- numeric(length(time))
dhealthy_dfc <- numeric(length(time))
dcancer_dfc <- numeric(length(time))

# Get initial locations
data <- subset(sim_data, V1 == 0.00)
healthy0 <- subset(data, V3 == "H")[,4:6]
cancer0 <- subset(data, V3 == "C")[,4:6]

i <- 1
for (t in time) {
  data <- subset(sim_data, V1 == t)
  healthy <- subset(data, V3 == "H")
  cancer <- subset(data, V3 == "C")
  healthy_dfc[i] <- mean(mapply(function(x, y, z){x^2+y^2+z^2},
                                healthy$V4, healthy$V5, healthy$V6))
  cancer_dfc[i] <- mean(mapply(function(x, y, z){x^2+y^2+z^2},
                               cancer$V4, cancer$V5, cancer$V6))
  i <- i + 1
}
plot(time, healthy_dfc, col="blue")
plot(time, cancer_dfc, col="red")
