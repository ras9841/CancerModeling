# File:   calc_dfc.r
# Author: Roland Sanford (ras9841@rit.edu)
# Description:
#   Plots the distance-from-center of the healthy (blue) and cancerous (red) cells.

library(data.table)
library(ggplot2)


# Read in data
base <- "~/Documents/CancerModeling/outputs/diff_surf-"
data_location <- paste(base, 1, "_msd.csv", sep="")
sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)
num_sim <- 2

time = sim_data$V1
healthy_msd = numeric(length=length(time))
cancer_msd = numeric(length=length(time))

for (i in 1:num_sim) {
  data_location <- paste(base, i, "_msd.csv", sep="")
  sim_data <- fread(data_location, sep="\t", header=FALSE, skip=4)
  healthy_msd <- healthy_msd + sim_data$V2
  cancer_msd <- cancer_msd + sim_data$V3
}

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