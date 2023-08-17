################# use this
library(dplyr)
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


###Import files in most base form to do only minor edits then feed in here - change only directory 
#and chamge title at the bottom 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/SFP10/")
combined=read.csv("Gordoniaall_confirmatory.csv")
#read in all columns
colnames(combined) <- c("X","variable","value","shape","shape2","tool","pvalue")
#but only select these, leaveoutwish plvaues intil upset plot
combined=combined[,c("X","variable","value","shape","shape2","tool")]


combined <- ddply(combined, .(X), transform, rescale = rescale(value))
combined <- ddply(combined, .(variable), transform, rescale2 = rescale(value))
combined <- ddply(combined, .(tool), transform, rescale3 = rescale(value))


combined$variable=sub("^(\\S+) (\\S+) ", "\\1 \\2)~(", combined$variable)
combined$variable=paste0("italic(", combined$variable)
combined$variable=sub(" ", "", combined$variable)
combined$variable <- paste0(combined$variable, ")")


#groups x axis by tool
combined$rename <- paste(combined$tool,combined$X)
#groups x axis by phage
#combined$rename <- paste(combined$X,combined$tool)


#combined <- combined[order(combined$X),]

labelsf=c("GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8")


library(forcats)
library(ggplot2)
# Plot.
q = ggplot() +
  geom_point(data=combined, aes(rename, fct_inorder(variable), color=rescale3, shape = shape), size = 5) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_x_discrete(labels = labelsf) + 
  scale_y_discrete(labels = parse(text = unique(combined$variable))) + 
  #keep sizes 18 and 14 for GMA/GRU graphs 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
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
  geom_point(data=combined, aes(rename, fct_inorder(variable), shape = shape2), size = 5, stroke=2) +
  
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_x_discrete(labels = labelsf) + 
  scale_y_discrete(labels = parse(text = unique(combined$variable))) +
  
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

q
#maybe put title on hold
#q + ggtitle("SFP10") +
#  theme(plot.title = element_text(size = 22, hjust = 0.5))


################################################################
#############for ecoli viruses
###########################################
###Import files in most base form to do only minor edits then feed in here - change only directory 
#and chamge title at the bottom 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/KFS-EC3/")
combined=read.csv("KFS-EC3_confirmatory.csv")

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/HY01/")
combined=read.csv("HY01_confirmatory.csv")

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/SFP10/")
combined=read.csv("SFP10_Onlyconfirmatory.csv")


#read in all columns
colnames(combined) <- c("X","variable","value","shape","shape2","tool","pvalue")
#but only select these, leaveoutwish plvaues intil upset plot
combined=combined[,c("X","variable","value","shape","shape2","tool")]


#combined <- ddply(combined, .(X), transform, rescale = rescale(value))
#combined <- ddply(combined, .(variable), transform, rescale2 = rescale(value))
#rescale by tool
combined <- ddply(combined, .(tool), transform, rescale3 = rescale(value))


combined$variable=sub("^(\\S+) (\\S+) ", "\\1 \\2)~(", combined$variable)
combined$variable=paste0("italic(", combined$variable)
combined$variable=sub(" ", "", combined$variable)
combined$variable <- paste0(combined$variable, ")")



library(forcats)
library(ggplot2)
# Plot.
q = ggplot() +
  geom_point(data=combined, aes(tool, fct_inorder(variable), color=rescale3, shape = shape), size = 12) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = unique(combined$variable))) + 
  #keep sizes 18 and 14 for GMA/GRU graphs 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 15), 
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
  geom_point(data=combined, aes(tool, fct_inorder(variable), shape = shape2), size = 11, stroke=2) +
  
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = unique(combined$variable))) +
  
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

q
#maybe put title on hold
#q + ggtitle("SFP10") +
#  theme(plot.title = element_text(size = 22, hjust = 0.5))

