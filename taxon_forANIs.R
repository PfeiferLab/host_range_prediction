setwd("/Users/pfeiferlab/Documents/hostrange/")
names=read.csv("nameslist.csv", header=FALSE)
setwd("/Users/pfeiferlab/")
library(taxonomizr)
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

#for gordonia phages host dataset
php2=php[grepl("Gordonia", php$V9),]
php_gonly=php2[grepl("malaquae|terrae|hydrophobica|rubripertincta", php2$V9),]

vhmn2=vhmn[grepl("Gordonia", vhmn$V9),]
vhmn_gonly=vhmn2[grepl("malaquae|terrae|hydrophobica|rubripertincta", vhmn2$V9),]

to=rbind(php_gonly, vhmn_gonly)
to2=to[,c("V1")]
setwd("/Users/pfeiferlab/Documents/hostrange/")
write.table(to2, "outfile.txt",row.names=FALSE,quote=FALSE,col.names=FALSE)

##use this outfile - no looping needed 
###bit-dl-ncbi-assemblies -w outfile.txt -f fasta

#all other e coli viruses 
#oops... younwant ll the ones it infects... not all the ones you have datasets for...
php3=php[grepl("Salmonella|Shigella|Escherichia", php$V9),]
php_sonly=php3[grepl("enterica|flexneri|coli|enterocolitica", php3$V9),]

vhmn3=vhmn[grepl("Salmonella|Shigella|Escherichia", vhmn$V9),]
vhmn_sonly=vhmn3[grepl("enterica|flexneri|coli|enterocolitica", vhmn3$V9),]

#SO - using this check 
#to=rbind(php_sonly, vhmn_sonly)
#to2=to[,c("V1")]
#str_detect(to$V1, to$V1) - it seems all genomes in the VHMN for the e coli hosts dataset (vhmn-sonly = 1002)
#are present in the php dataset (8960) - so i see no need to include other entires  - let use only the php dataset 
php_sonly1=php_sonly[,c("V1")]

vhmn_sonly1=vhmn_sonly[,c("V1")]

#PHP has 8000 files, use the VHMN dataset and explain caveat 

#str_detect(to$V1, to$V1)
setwd("/Users/pfeiferlab/Documents/hostrange/")
write.table(vhmn_sonly1, "outfile_s.txt",row.names=FALSE,quote=FALSE,col.names=FALSE)

####once youve run anvio plot and fix the labels
library(ggplot2)
library(reshape2)
library(viridis)

setwd("/Users/pfeiferlab/Documents/hostrange/ani_gordoniadb_long_output4/")
ani.dat <- read.csv("ANIb_percentage_identity.txt", sep="\t", header = T)


#if your entries include the .X part need to remove version=base
#taxaId<-taxonomizr::accessionToTaxa(hosts_list,"accessionTaxa.sql", version='base')
setwd("/Users/pfeiferlab/")
taxaId<-taxonomizr::accessionToTaxa(ani.dat$key,"accessionTaxa.sql")
output=as.data.frame(getTaxonomy(taxaId,'accessionTaxa.sql'))

#ok this doesnt work (as well) because it also detects a duplicate anytime "gordonia" pops up multiple times
#but it does quick chekc that the first 8 are duplicates so I'm just going to amnually scrub those 
#library(stringr)
#str_detect(ani.dat$key, ani.dat$key)


library("dplyr")
ani.dat$rename=ifelse(is.na(output$species), ani.dat$key, paste(output$species,ani.dat$key))
ani.dat=select(ani.dat, -c(key)) 
colnames(ani.dat)=ani.dat$rename

#drop na column at the end or you cant use the sort/reorder functions 
ani.dat<- ani.dat [1: ncol(ani.dat)-1 ]
#drop first 7 columns because theyre repeats 
df <- ani.dat[ -c(1:7) ]

ani.dat=ani.dat[,order(colnames(ani.dat))]

#reorder so the tested genomes are in the bottom left corner 
ani.dat <- ani.dat %>%
  select(Gordonia_hydrophobica_44015, Gordonia_malaquae_44454,
         Gordonia_malaquae_44464, Gordonia_rubripertincta_DSM_43197, Gordonia_terrae_43249, everything())

#select only first 5
#this is speciically for the evolution talk nice zoomed in cersion
#ani.dat=ani.dat[1:5]

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
  
  scale_y_discrete(labels=c("G. hydrophobica (44015)","G. malaquae (44454)","G. malaquae (44464)",
                                             "G. rubripertinca (43197)","G. terrea (43249)")) +
  scale_x_discrete(labels=c("G. hydrophobica (44015)","G. malaquae (44454)","G. malaquae (44464)",
                            "G. rubripertinca (43197)","G. terrea (43249)")) 


#same thing for ecoli but setup is alightly different due to how the ATCC contigs are labelled 
setwd("/Users/pfeiferlab/Documents/hostrange/ani_ecoli_output2/")
ani.dat <- read.csv("ANIb_percentage_identity.txt", sep="\t", header = T)
#drop key column - its different from the manual columns i just fixed 
ani.dat=select(ani.dat, -c(key)) 
ani.dat$key=colnames(ani.dat)

setwd("/Users/pfeiferlab/")
library(taxonomizr)
taxaId<-taxonomizr::accessionToTaxa(ani.dat$key,"accessionTaxa.sql")
output=as.data.frame(getTaxonomy(taxaId,'accessionTaxa.sql'))

ani.dat$rename=ifelse(is.na(output$species), ani.dat$key, paste(output$species,ani.dat$key))

ani.dat=select(ani.dat, -c(key)) 
colnames(ani.dat)=ani.dat$rename
#drop weird hanging NA column leftover from renaming
ani.dat<- ani.dat [1: ncol(ani.dat)-1 ]

#reodrer all alphabetically first 
ani.dat=ani.dat[,order(colnames(ani.dat))]

colnames(ani.dat)

#then reorder so the tested genomes are in the bottom left corner 
ani.dat <- ani.dat %>%
  select(Escherichia_coli_ATCC_43888, Escherichia_coli_ATCC_15144,	Escherichia_coli_ATCC_10536, 
         Escherichia_coli_ATCC_43890, Escherichia_coli_ATCC_BAA_2196,	Escherichia_coli_ATCC_BAA_2192,
         Escherichia_coli_ATCC_43894,	Escherichia_coli_ATCC_25922,	Escherichia_coli_ATCC_43895, Escherichia_coli_ATCC_35150,	
         Escherichia_colistr.K12.MG1655_source29,
         Shigella_flexneri_ATCC_29903, Shigella_sonnei_ATCC_9290,	Shigella_flexneri2astr.2457T,Shigella_flexneri_ATCC_12022,
         Salmonella_enterica_subsp_enterica_ATCC_13076,		
         Salmonella_enterica_subsp_enterica_ATCC_13311,
         Salmonella_enterica_ATCC_14028,	
         Salmonella_entericasubsp_Typhimuriumstr.SL1344,
         Salmonella_entericaserovar_Typhimuriumstr.LT2, everything())


data = cor(ani.dat[sapply(ani.dat, is.numeric)])
data1 <- melt(data)


#this is NOT the lenght of data4 but the lenght of ANI.dat
breaks=levels(data1$Var1)[c(1,12,16,21,114,197)]

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
  
  
  
 scale_y_discrete(breaks = breaks, labels=c("E.coli ATCC","Shigella ATCC","Salmonella ATCC",
                                            "E.coli database","Shigella database","Salmonella database")) +
  scale_x_discrete(breaks = breaks, labels=c("E.coli ATCC","Shigella ATCC","Salmonella ATCC",
                                             "E.coli database","Shigella database","Salmonella database")) 


#ok now just the confirmatory guys zoomed in
#then reorder so the tested genomes are in the bottom left corner 
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
