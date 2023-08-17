
# Load libraries.
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(scales)
library(dplyr)
library(RColorBrewer)
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/")

####correlation plot
comb=read.csv("hy01_phirbo_gordonia_paper.list.matrix.melted.contignumber.csv")

par(mfrow=c(2,2))  


# Creating the plot
plot(comb$contignum, comb$phirbovalue, pch = 19, col = "lightblue",
     ylab="Phirbo Score",
     xlab="Contig Number")
# Regression line
abline(lm(comb$phirbovalue ~ comb$contignum), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(comb$contignum, comb$phirbovalue), 2)), x = 80, y = 0.05)


# Creating the plot
plot(comb$contignum, comb$wishvalue, pch = 19, col = "lightblue",
     ylab="WIsH Score",
     xlab="Contig Number")
# Regression line
abline(lm(comb$wishvalue ~ comb$contignum), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(comb$contignum, comb$wishvalue), 2)), x = 80, y = 0.05)

# Creating the plot
plot(comb$contignum, comb$phpvalue, pch = 19, col = "lightblue",
     ylab="PHP Score",
     xlab="Contig Number")
# Regression line
abline(lm(comb$phpvalue ~ comb$contignum), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(comb$contignum, comb$phpvalue), 2)), x = 80, y = 1450)


# Creating the plot
plot(comb$contignum, comb$vhmvalue, pch = 19, col = "lightblue",
     ylab="VHM Score",
     xlab="Contig Number")
# Regression line
abline(lm(comb$vhmvalue ~ comb$contignum), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(comb$contignum, comb$vhmvalue), 2)), x = 80, y = 0.05)













