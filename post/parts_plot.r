# File:   parts_plot.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the positions of the healthy (blue) and cancerous (red) cells within 
#   the bounding sphere.

library(rgl)
N <- 256
Rb <-0.0142276  # taken from debug output
base <- "~/Documents/CancerModeling/outputs/long_test-"
data_location <- paste(base, 1, ".csv", sep="")
sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)

start = dim(sim_data)[1] - N
while (start < N) {
    data <- sim_data[start:(start+N-1),]
    healthy <- subset(data, V3 == "H")
    cancer <- subset(data, V3 == "C")
    plot3d(healthy$V4, healthy$V5, healthy$V6, xlab="", ylab="", zlab="",
        col="blue", size=2, type="s", axes=FALSE, add=FALSE)
    plot3d(cancer$V4, cancer$V5, cancer$V6, col="red", size=2, type="s", add=TRUE)
    rgl.spheres(x=0, y=0, z=0, radius = Rb, alpha=.3, add=TRUE)
    start <- start+N
    Sys.sleep(2)
}