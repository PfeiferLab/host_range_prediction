## Written by: Abby Howell
## Modified by: Cyril Versoza
## Date: August 28, 2023

# Load necessary libraries:
library(taxonomizr)
library(stringr)
library(ggplot2)
library(reshape2)
library(viridis)
library(taxonomizr)
library("dplyr")
library(gtools)

# Use nameslist.csv to rename downloaded accession numbers:
names=read.csv("nameslist.csv", header=FALSE)
taxaId<-taxonomizr::accessionToTaxa(names$V1,"accessionTaxa.sql")
output=as.data.frame(getTaxonomy(taxaId,'accessionTaxa.sql'))

output$name=names$V1
output2=output[grepl("Escherichia|Entero|Gordonia|Nocardia|Rhodococcus|Salmonella|
                     Shigella|Aeromonas|Bac", output$species),]
setwd("/Users/pfeiferlab/Documents/hostrange/")
write.csv(output2, "remove.names.csv")

setwd("/Users/pfeiferlab/Downloads/")
php=read.csv("60105Taxonomy.txt",header=FALSE, sep="\t")
vhmn=read.csv("hostTaxa_VHMN.csv",header=FALSE)

# For Gordonia phages host dataset
php2=php[grepl("Gordonia", php$V9),]
php_gonly=php2[grepl("malaquae|terrae|hydrophobica|rubripertincta", php2$V9),]

vhmn2=vhmn[grepl("Gordonia", vhmn$V9),]
vhmn_gonly=vhmn2[grepl("malaquae|terrae|hydrophobica|rubripertincta", vhmn2$V9),]

to=rbind(php_gonly, vhmn_gonly)
to2=to[,c("V1")]
write.table(to2, "outfile.txt",row.names=FALSE,quote=FALSE,col.names=FALSE)

## Use this outfile - no looping needed 
### bit-dl-ncbi-assemblies -w outfile.txt -f fasta

# All other E. coli viruses 
php3=php[grepl("Salmonella|Shigella|Escherichia", php$V9),]
php_sonly=php3[grepl("enterica|flexneri|coli|enterocolitica", php3$V9),]

vhmn3=vhmn[grepl("Salmonella|Shigella|Escherichia", vhmn$V9),]
vhmn_sonly=vhmn3[grepl("enterica|flexneri|coli|enterocolitica", vhmn3$V9),]

php_sonly1=php_sonly[,c("V1")]

vhmn_sonly1=vhmn_sonly[,c("V1")]

write.table(vhmn_sonly1, "outfile_s.txt",row.names=FALSE,quote=FALSE,col.names=FALSE)

## Average nucleotide identity:
ani.dat <- read.csv("ANIb_percentage_identity.txt", sep="\t", header = T)

# If your entries include the .X part need to remove version=base
# taxaId<-taxonomizr::accessionToTaxa(hosts_list,"accessionTaxa.sql", version='base')
taxaId<-taxonomizr::accessionToTaxa(ani.dat$key,"accessionTaxa.sql")
output=as.data.frame(getTaxonomy(taxaId,'accessionTaxa.sql'))

ani.dat$rename=ifelse(is.na(output$species), ani.dat$key, paste(output$species,ani.dat$key))
ani.dat=select(ani.dat, -c(key)) 
colnames(ani.dat)=ani.dat$rename

# Drop na column at the end or you cant use the sort/reorder functions 
ani.dat<- ani.dat [1: ncol(ani.dat)-1 ]
# Drop first 7 columns because theyre repeats 
ani.dat <- ani.dat[ -c(1:7) ]

ani.dat=ani.dat[,order(colnames(ani.dat))]
ani.dat=ani.dat[,mixedsort(colnames(ani.dat))]

# Reorder so the tested genomes are in the bottom left corner 
ani.dat <- ani.dat %>%
  select(Gordonia_malaquae_44454,
         Gordonia_malaquae_44464, Gordonia_rubripertincta_DSM_43197, Gordonia_terrae_43249, Gordonia_hydrophobica_44015,
         everything())

data = cor(ani.dat[sapply(ani.dat, is.numeric)])
data1 <- melt(data)

# Create heatmap
ggplot(data1, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_viridis() + # nice color scheme!ß
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), 
        axis.text.y = element_text(size=10), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  labs(fill="ANI") + 
# Don't plot the discrete axis labels part unless youre only doing the zooomed in 5 
  scale_y_discrete(labels=c("G. hydrophobica (44015)","G. malaquae (44454)","G. malaquae (44464)",
                                             "G. rubripertinca (43197)","G. terrea (43249)")) +
  scale_x_discrete(labels=c("G. hydrophobica (44015)","G. malaquae (44454)","G. malaquae (44464)",
                            "G. rubripertinca (43197)","G. terrea (43249)")) 


# Same thing for ecoli but setup is alightly different due to how the ATCC contigs are labelled 

ani.dat <- read.csv("ANIb_percentage_identity.txt", sep="\t", header = T)
# Drop key column - its different from the manual columns i just fixed 
ani.dat=select(ani.dat, -c(key)) 
ani.dat$key=colnames(ani.dat)

taxaId<-taxonomizr::accessionToTaxa(ani.dat$key,"accessionTaxa.sql")
output=as.data.frame(getTaxonomy(taxaId,'accessionTaxa.sql'))

ani.dat$rename=ifelse(is.na(output$species), ani.dat$key, paste(output$species,ani.dat$key))

ani.dat=select(ani.dat, -c(key)) 
colnames(ani.dat)=ani.dat$rename
# Drop hanging NA column leftover from renaming
ani.dat<- ani.dat [1: ncol(ani.dat)-1 ]

# Reorder all alphabetically first.
ani.dat=ani.dat[,order(colnames(ani.dat))]

colnames(ani.dat)
ani.dat=ani.dat[,mixedsort(colnames(ani.dat))]

# Then reorder so the tested genomes are in the bottom left corner.
ani.dat <- ani.dat %>%
  select(Escherichia_coli_ATCC_10536,
        Escherichia_coli_ATCC_35150,
        Escherichia_coli_ATCC_43888, 
        Escherichia_coli_ATCC_43890,
        Escherichia_coli_ATCC_43894,
        Escherichia_coli_ATCC_43895,
        Salmonella_enterica_subsp_enterica_ATCC_13076,		
        Salmonella_enterica_ATCC_14028,	
        Salmonella_entericaserovar_Typhimuriumstr.LT2,
        Salmonella_entericasubsp_Typhimuriumstr.SL1344,
        Shigella_flexneri_ATCC_12022,
        Shigella_flexneri_ATCC_29903,
        Shigella_flexneri2astr.2457T,
        Shigella_sonnei_ATCC_9290,
        everything())


data = cor(ani.dat[sapply(ani.dat, is.numeric)])
data1 <- melt(data)

ani.dat = ani.dat[, -which(names(ani.dat) %in% c("Escherichia_coli_ATCC_15144", "Escherichia_coli_ATCC_25922",
                                    "Escherichia_coli_ATCC_BAA_2192","Escherichia_coli_ATCC_BAA_2196",
                                    "Escherichia_colistr.K12.MG1655_source29",
                                    "Salmonella_enterica_subsp_enterica_ATCC_13311",
                                    "Shigella flexneri NC_008258.1"))]

# Show first instance of Escherichia, Salmonella, and Shigella.
breaks=levels(data1$Var1)[c(1,7,11,15,108,191)]


# Create heatmap
ggplot(data1, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_viridis() + # nice color scheme!ß
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), # make the x-axis labels vertical for readability
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  labs(fill="ANI") + 
  scale_y_discrete(breaks = breaks) +
  scale_x_discrete(breaks = breaks) 

# Ok now just the confirmatory guys zoomed in
# Then reorder so the tested genomes are in the bottom left corner 
ani.dat2 = ani.dat[1:20]

data = cor(ani.dat2[sapply(ani.dat2, is.numeric)])
data1 <- melt(data)

# Create heatmap
ggplot(data1, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_viridis() + # nice color scheme!ß
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), # make the x-axis labels vertical for readability
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  labs(fill="ANI") 