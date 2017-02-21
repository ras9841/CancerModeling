# File:   calc_dfc.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the distance-from-center of the healthy (blue) and cancerous (red) cells.

library(doParallel)
library(foreach)
library(data.table)
library(ggplot2)

calc_sdfc <- function (xvals, yvals, zvals) {
  vals = numeric(length(xvals))
  for (i in 1:length(xvals))  {
    vals[i] <- xvals[i]^2 + yvals[i]^2 + zvals[i]^2
  }
  vals
}

get_healthy_msdfc <- function(time, sim_data){
  foreach(i = 1:length(time), .combine = c, .export = "calc_sdfc")  %dopar%  {
    start <- 1+256*(i-1)
    end <- start + 256
    h_index = start:(start+128-1)
    healthy <- sim_data[h_index, 4:6];
    mean(calc_sdfc(healthy[,1], healthy[,2], healthy[,3]))
  }
}

get_cancer_msdfc <- function(time, sim_data){
  foreach(i = 1:length(time), .combine = c, .export = "calc_sdfc")  %dopar%  {
    start <- 1+256*(i-1)
    end <- start + 256
    c_index = (start+128):end;
    cancer <- sim_data[c_index, 4:6];
    mean(calc_sdfc(cancer[,1], cancer[,2], cancer[,3]))
  }
}

main <- function() {
  # Read in data
  data_location <- "~/Documents/CancerModeling/outputs/testC-1.csv"
  sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)

  # Setup time and MSD vectors
  time <- unique(sim_data$V1)

  # Set up threads
  no_cores <- 4;
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  # calculate MSDFC
  healthy_dfc <- get_healthy_msdfc(time, sim_data)
  cancer_dfc <- get_cancer_msdfc(time, sim_data)
  stopCluster(cl)

  # Plot Results
  df <- data.frame(time,healthy_dfc,cancer_dfc)
  g <-ggplot(df, aes(time)) +                  
    geom_line(aes(y=healthy_dfc), colour="blue") + 
    geom_line(aes(y=cancer_dfc), colour="red") +
    ylab("MSDFC (sigma^2)") +
    xlab("Time (tau)")
  g
}

main()