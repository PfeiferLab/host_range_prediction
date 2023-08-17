#same thing for ecoli but setup is alightly different due to how the ATCC contigs are labelled 
#actually here on lab desktop 
#setwd("/Users/pfeiferlab/Documents/hostrange/ani_ecoli_hosts_seq_output/")
#but temp working on it locally 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/")
library(stringr)
library("dplyr")

ani.dat <- read.csv("ANIb_percentage_identity.txt", sep="\t", header = T)
#pain in the ass renaming schemme 
#get into correct format first from echo $x and head-n1

key=read.csv("ecoli_hosts_name_key.csv",header=FALSE)
key2=as.data.frame(matrix(key$V1, ncol = 2, byrow = TRUE))

key2$V1=sub('.fasta.long', '', key2$V1)
key2$V1=sub('.fna.long', '', key2$V1)

key2$V2=sub('>', '', key2$V2)

#get rid of stupid X in column names
colnames(ani.dat)=sub('X', '', colnames(ani.dat))


ani.dat$key <- dplyr::recode(
  ani.dat$key, 
  !!!setNames(as.character(key2$V1), as.character(key2$V2))
)

colnames(ani.dat) <- dplyr::recode(
  colnames(ani.dat), 
  !!!setNames(as.character(key2$V1), as.character(key2$V2))
)

####
####DO NOT USE THIS UNLESSS SUBSETTING FOR A REASON
#why does salmonella v ecoli look less different here than the alone version?
#drop rows first then cols
php3=ani.dat[grepl("Salmonella|Shigella|Escherichia", ani.dat$key),]
php4=php3[,grepl("Salmonella|Shigella|Escherichia", colnames(php3))]
#order alphabeticLLY
ani.dat=php4[,order(colnames(php4))]
###################
#ok guess it really does matter whos being comapred to what its all on a relative scale 


###STRAIGHT TO THIS FOR CONFIRMATORY ECOLI DATASET
#order alphabeticLLY
ani.dat=ani.dat[,order(colnames(ani.dat))]

library(gtools)
ani.dat=ani.dat[,mixedsort(colnames(ani.dat))]


#order so all actual hosts are in the corner
ani.dat <- ani.dat %>%
  select(Escherichia_coli_ATCC_10536, Escherichia_coli_ATCC_35150, 
         Escherichia_coli_ATCC_43888, Escherichia_coli_ATCC_43890, Escherichia_coli_ATCC_43894,
         Escherichia_coli_ATCC_43895, Escherichiacolistr.K12substr.MG1655,
        Salmonella_enterica_ATCC_13076,
        Salmonella_Typhimuriumstr.LT2,
         Salmonella_Typhimuriumstr.SL1344,
         Salmonella_enterica_ATCC_13076,		
         Salmonella_enterica_ATCC_14028,
         Shigella_flexneri_ATCC_29903,
        Shigella_flexneri2astr.2457T,
        Shigella_sonnei_ATCC_9290, everything())




library(reshape2)
library(ggplot2)
library(viridis)

data = cor(ani.dat[sapply(ani.dat, is.numeric)])
data1 <- melt(data)

# Create heatmap
ggplot(data1, aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) +
  scale_fill_viridis() + # nice color scheme!ß
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  labs(fill="ANI") 


