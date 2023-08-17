# Written by: Cyril Joseph Versoza
# Date: March 26, 2022

# This script will be used to generate a heatmap from the ANI values calculated by DNA Master.
# Specifically, ANI values between BiggityBass and other Gordonia phages.

# Load libraries.
library(ggplot2)
library(reshape)
library(viridis)

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/")
# Load table.
#ani.dat <- read.csv("ClusterPs_averagenucleotidedensity.csv", header = T)
ani.dat <- read.csv("ANIb_percentage_identity.names.txt", sep="\t", header = T)

names=read.csv("outfile_names.csv",header=FALSE)
colnames(ani.dat) = names$V9
names=read.csv("outfile_names.csv",header=TRUE)
ani.dat$V9=names$V9

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/")
ani.dat <- read.csv("ANIb_percentage_identity_viruses.txt", sep="\t", header = T)
names=read.csv("ani_nameslist.csv",header=FALSE)

ani.dat=select(ani.dat, -c(key)) 
colnames(ani.dat) = names$V1
ani.dat$key=colnames(ani.dat)



data = cor(ani.dat[sapply(ani.dat, is.numeric)])

data1 <- melt(data)

# Create heatmap
ggplot(data1, aes(X1, X2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_viridis() + # nice color scheme!
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  labs(fill="ANI") 


library("d3heatmap")
d3heatmap(data, colors = "RdYlBu")

library("d3heatmap")
d3heatmap(scale(mtcars), colors = "RdYlBu",Rowv=FALSE,Colv = FALSE)

library(taxonomizr)
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange")
hosts=read.csv("hmp_host_list.txt",header=FALSE)
prepareDatabase('accessionTaxa.sql',resume=FALSE)

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/SEA-Phages/")

taxaId<-taxonomizr::accessionToTaxa(hosts,"accessionTaxa.sql")


output=getTaxonomy(taxaId,'accessionTaxa.sql')


