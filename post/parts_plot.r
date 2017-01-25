# File:   parts_plot.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the positions of the healthy (blue) and cancerous (red) cells within 
#   the bounding sphere.

library(rgl)
Rb <-0.0142276  # taken from debug output
data_location <- "~/Documents/CancerModeling/outputs/testA.csv"
sim_data <- read.csv(data_location, comment.char="#", sep="\t", header=FALSE)
healthy <- subset(sim_data, V3 == "H")
cancer <- subset(sim_data, V3 == "C")

plot3d(healthy$V4, healthy$V5, healthy$V6, xlab="", ylab="", zlab="",
       col="blue", size=2, type="s", axes=FALSE)
plot3d(cancer$V4, cancer$V5, cancer$V6, col="red", size=2, type="s", add=TRUE)
rgl.spheres(x=0, y=0, z=0, radius = Rb, alpha=.3, add=TRUE)
