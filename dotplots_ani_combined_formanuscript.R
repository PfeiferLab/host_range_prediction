library(dplyr)
library(data.table)
library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(scales)
library(dplyr)
library(RColorBrewer)

#all gordonia and e. coli phage datasets combined
combined=read.csv("everything_all.csv")

#read in all columns
colnames(combined) <- c("X","variable","value","shape","shape2","tool","pvalue")
#but only select these, leaveoutwish plvaues intil upset plot
#combined=combined[,c("X","variable","value","shape","shape2","tool")]
#combined <- ddply(combined, .(X), transform, rescale = rescale(value))
#combined <- ddply(combined, .(variable), transform, rescale2 = rescale(value))

combined$finalvalue=ifelse(combined$tool=="WIsH", combined$pvalue, combined$value)
combined$finalvalue2 = ifelse(combined$tool=="VHM" | combined$tool=="WIsH", 1-combined$finalvalue,combined$finalvalue)

#rescale per tool across all datasets
combined <- ddply(combined, .(tool), transform, rescale3 = rescale(finalvalue2))

combined$variable=sub("^(\\S+) (\\S+) ", "\\1 \\2)~(", combined$variable)
combined$variable=paste0("italic(", combined$variable)
combined$variable=sub(" ", "", combined$variable)
combined$variable <- paste0(combined$variable, ")")

table(combined$tool)

#groups x axis by tool
combined$rename <- paste(combined$tool,combined$X)

GMAs=combined[combined$X %like% "G",]
labelsf=c("GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8",
          "GMA2","GMA3","GMA4","GMA5","GMA6","GMA7","GRU1","GRU3","GTE2","GTE5","GTE6","GTE7","GTE8")


library(forcats)
library(ggplot2)
# Plot.
q = ggplot() +
  geom_point(data=GMAs, aes(rename, fct_inorder(variable), color=rescale3, shape = shape), size = 5) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_x_discrete(labels = labelsf) + 
  scale_y_discrete(labels = parse(text = unique(GMAs$variable))) + 
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
  
  ###this whole next part is just to get the border
  geom_point(data=GMAs, aes(rename, fct_inorder(variable), shape = shape2), size = 5, stroke=2) +
  
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_x_discrete(labels = labelsf) + 
  scale_y_discrete(labels = parse(text = unique(GMAs$variable))) +
  
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
ecoli=combined[combined$X == "KFS-EC3",]
ecoli=combined[combined$X == "HY01",]
ecoli=combined[combined$X == "SFP10",]


library(forcats)
library(ggplot2)
# Plot.
q = ggplot() +
  geom_point(data=ecoli, aes(tool, fct_inorder(variable), color=rescale3, shape = shape), size = 12) +
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = unique(ecoli$variable))) + 
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
  
  ###this whole next part is just to get the border
  geom_point(data=ecoli, aes(tool, fct_inorder(variable), shape = shape2), size = 11, stroke=2) +
  
  #scale_shape_identity() +
  #scale_color_viridis(option = "viridis", name = "log likelihood") +
  scale_y_discrete(labels = parse(text = unique(ecoli$variable))) +
  
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

