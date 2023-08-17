# Written by: Cyril Joseph Versoza
# Date: May 6, 2022

# This script is for plotting the Phirbo results with the added variable to identify which host was experimentall validated.

# Load libraries.
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(scales)
library(dplyr)
library(RColorBrewer)

# Load new dataframe.
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/")

#this one you dont need to write out and redo just write shape column to it directly
paper.wish.dat <- read.csv("kfs-ec3_wish_prediction.list",sep="\t")
paper.wish.dat <- read.csv("hy01_wish_prediction.list",sep="\t")
paper.wish.dat <- read.csv("sfp10_wish_prediction.list",sep="\t")


colnames(paper.wish.dat) <- c("X","variable","LogLikelihood","pvalue","shape")

#paper.wish.dat <- read.csv("wish_gordonia_paper.shape.csv",sep=",")
#colnames(paper.wish.dat) <- c("X","variable","LogLikelihood","pvalue","shape")
# Rescale values.
paper.wish.shape.dat <- ddply(paper.wish.dat, .(X), transform, rescale = rescale(LogLikelihood))
paper.wish.shape.dat <- ddply(paper.wish.dat, .(variable), transform, rescale2 = rescale(LogLikelihood))

#add phirbo col
paper.melt.shape.dat.phirbo <- read.csv("sfp10_phirbo_gordonia_paper.list.matrix.melted.csv")
#add virhostmatcher col
paper.vhm.shape.dat.vir <- read.csv("virhostmatcher_gordonia_paper.melted.csv",sep=",")
colnames(paper.vhm.shape.dat.vir) <- c("X","variable","value","shape")

wish=paper.wish.shape.dat[,c("X","variable","LogLikelihood","shape")]
colnames(wish) <- c("X","variable","value","shape")

phirbo=paper.melt.shape.dat.phirbo[,c("X","variable","value","shape")]

df1=rbind(wish,phirbo)
df2=rbind(df1,paper.vhm.shape.dat.vir)


library(tidyverse)
php=read.csv("HY01_PHPHostKmer_Prediction_Allhost.csv",sep=",")
php <- melt(php)
write.csv(php,"HY01_PHP.melted.shape.csv", row.names = FALSE)

php=read.csv("SFP10_PHPHostKmer_Prediction_Allhost.csv",sep=",")
php <- melt(php)
write.csv(php,"SFP10_PHP.melted.shape.csv", row.names = FALSE)


php=read.csv("KFS-EC3_PHPHostKmer_Prediction_Allhost.csv",sep=",")
php <- melt(php)
write.csv(php,"KFS-EC3_PHP.melted.shape.csv", row.names = FALSE)





setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/")

####correlation plot??????
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















VHM=read.csv("hy01_virhostmatcher_gordonia_paper.melted.contignumber.csv")
VHM$N50=VHM$lengthlongestcontig / VHM$total.length
plot(VHM$N50, VHM$value)

wish=read.csv("hy01_wish_prediction_contignumb.csv")
plot(wish$contignum, wish$value)

################# use this
library(dplyr)

par(mfrow=c(3,1))  


###KFS-EC3
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH")
virhostmatcher=read.csv("KFS-EC3_virhostmatcher_gordonia_paper.melted.csv")
virhostmatcher$tool="VHM"
colnames(virhostmatcher) <- c("X","variable","value","shape","tool")
phirbo=read.csv("KFS-EC3_phirbo_gordonia_paper.list.matrix.melted.csv")
phirbo$tool="Phirbo"
phirbo=phirbo[,c("X","variable","value","shape","tool")]
wish=read.csv("kfs-ec3_wish_prediction.csv")
wish$tool="WIsH"
wish=wish[,c("Phage","Best.hit.among.provided.hosts","LogLikelihood","shape","tool")]
colnames(wish) <- c("X","variable","value","shape","tool")
php=read.csv("KFS-EC3_PHP.melted.shape.csv",sep=",")
php$tool="PHP"
colnames(php) <- c("X","variable","value","shape","tool")


#HY01
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH")
virhostmatcher=read.csv("HY01_virhostmatcher_gordonia_paper.melted.csv")
virhostmatcher$tool="virhostmatcher"
colnames(virhostmatcher) <- c("X","variable","value","shape","tool")
phirbo=read.csv("HY01_phirbo_gordonia_paper.list.matrix.melted.csv")
phirbo$tool="phirbo"
phirbo=phirbo[,c("X","variable","value","shape","tool")]
wish=read.csv("hy01_wish_prediction.csv")
wish$tool="wish"
wish=wish[,c("?..Phage","Best.hit.among.provided.hosts","LogLikelihood","shape","tool")]
colnames(wish) <- c("X","variable","value","shape","tool")
php=read.csv("HY01_PHP.melted.shape.csv",sep=",")
php$tool="php"
colnames(php) <- c("X","variable","value","shape","tool")


#SPF10
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH")
virhostmatcher=read.csv("SFP10_virhostmatcher_gordonia_paper.melted.csv")
virhostmatcher$tool="virhostmatcher"
colnames(virhostmatcher) <- c("X","variable","value","shape","tool")
phirbo=read.csv("SFP10_phirbo_gordonia_paper.list.matrix.melted.csv")
phirbo$tool="phirbo"
phirbo=phirbo[,c("X","variable","value","shape","tool")]
wish=read.csv("SFP10_wish_prediction.csv")
wish$tool="wish"
wish=wish[,c("?..Phage","Best.hit.among.provided.hosts","LogLikelihood","shape","tool")]
colnames(wish) <- c("X","variable","value","shape","tool")
php=read.csv("SFP10_PHP.melted.shape.csv",sep=",")
php$tool="php"
colnames(php) <- c("X","variable","value","shape","tool")

###susanne wants you to try and combine and rescale
library(data.table)
# Load libraries.
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(scales)
library(dplyr)
library(RColorBrewer)

datlist=list(wish,php,virhostmatcher,phirbo)

combined=rbindlist(datlist)


combined <- ddply(combined, .(X), transform, rescale = rescale(value))
combined <- ddply(combined, .(variable), transform, rescale2 = rescale(value))
combined <- ddply(combined, .(tool), transform, rescale3 = rescale(value))


combined$variable=sub("^(\\S+) (\\S+) ", "\\1 \\2)~", combined$variable)
combined$variable=paste0("italic(", combined$variable)
combined$variable=sub(" ", "", combined$variable)

#sodumb have to write it out - no doing an if else column does not work
#OHMYGOD THE THING YOU HAD TO CHANGE WAS 15-22 AND 16-21 FOR THE FILL SHAPES
#WHY DID YOU NOT WRITE THAT DOWN
#write.csv(combined,"combined.csv")
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH")
combined2=read.csv("combined.csv")

#write.csv(combined,"combined_hy01.csv")
combined2=read.csv("combined_hy01.csv")

#write.csv(combined,"combined_sfp10.csv")
combined2=read.csv("combined_sfp10.csv")

library(ggplot2)
# Plot.
ggplot() +
  geom_point(data=combined, aes(tool, variable, color=rescale3, shape = shape), size = 13) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = combined2$variable)) + 

  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  scale_color_distiller(type = "seq",
                        palette = "Greys",
                        breaks=c(0,0.5,1),
                        labels=c("low","medium","high"),
                                 name = "host likelihood") + 
  scale_continuous_identity(aesthetics = 'shape',
                            guide = 'legend',
                            breaks = c(15.75,16.00),
                            labels = c("no","yes"),
                            name = "experimentally validated host?") + 

  ###this whole next part is just to get the stupid border
  geom_point(data=combined2, aes(tool, variable, shape = shape), size = 12, stroke=2) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = combined2$variable)) + 
  
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, hjust = 1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +

  scale_continuous_identity(aesthetics = 'shape',
                            guide = 'legend',
                            breaks = c(15.75,16.00),
                            labels = c("no","yes"),
                            name = "experimentally validated host?") 
