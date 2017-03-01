# File:   calc_dfc.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the distance-from-center of the healthy (blue) and cancerous (red) cells.

library(doParallel)
library(foreach)
library(data.table)
library(ggplot2)

calc_sd <- function (xv, yv, zv, x0, y0, z0) {
  vals = numeric(length(xv))
  for (i in 1:length(xv))  {
    vals[i] <- (xv[i]-x0[i])^2 + (yv[i]-y0[i])^2 + (zv[i]-z0[i])^2
  }
  vals
}

get_healthy_msd <- function(time, sim_data, healthy0){
  foreach(i = 1:length(time), .combine = c, .export = "calc_sd")  %dopar%  {
    start <- 1+256*(i-1)
    end <- start + 256
    h_index = start:(start+128-1)
    healthy <- sim_data[h_index, 4:6]
    mean(calc_sd(healthy[,1], healthy[,2], healthy[,3], 
                   healthy0[,1], healthy0[,2], healthy0[,3]))
  }
}

get_cancer_msd <- function(time, sim_data, healthy0){
  foreach(i = 1:length(time), .combine = c, .export = "calc_sd")  %dopar%  {
    start <- 1+256*(i-1)
    end <- start + 256
    h_index = (start+128):end
    healthy <- sim_data[h_index, 4:6]
    mean(calc_sd(healthy[,1], healthy[,2], healthy[,3], 
                 healthy0[,1], healthy0[,2], healthy0[,3]))
  }
}

# Read in data
base <- "~/Documents/CancerModeling/outputs/step_test-"
data_location <- paste(base, 1, ".csv", sep="")
sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)
num_sim <- 2
  
# Setup time and MSD vectors
time <- unique(sim_data$V1)
  
# Get initial locations
data <- subset(sim_data, V1 == time[1])
healthy0 <- subset(data, V3 == "H")[,4:6]
cancer0 <- subset(data, V3 == "C")[,4:6]
  
# Set up threads
no_cores <- 4;
cl<-makeCluster(no_cores)
registerDoParallel(cl)
healthy_msd = numeric(length(time))
cancer_msd = numeric(length(time))
for (i in 1:num_sim) {
  data_location <- paste(base, i, ".csv", sep="")
  sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)
  # calculate MSD
  healthy_msd <- healthy_msd + get_healthy_msd(time, sim_data, healthy0)
  cancer_msd <- cancer_msd + get_healthy_msd(time, sim_data, cancer0)
}
stopCluster(cl)
healthy_msd <- healthy_msd/num_sim
cancer_msd <- cancer_msd/num_sim

# Plot Results
df <- data.frame(time,healthy_msd,cancer_msd)
g <-ggplot(df, aes(time)) +                  
  geom_line(aes(y=healthy_msd), colour="blue") + 
  geom_line(aes(y=cancer_msd), colour="red") +
  ylab("MSD (sigma^2)") +
  xlab("Time (tau)")
g