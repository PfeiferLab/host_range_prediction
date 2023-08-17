setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/")

virhostmatcher=read.csv("virhostmatcher_gordonia_paper.melted.csv")
phirbo=read.csv("phirbo_gordonia_paper.list.matrix.melted.csv")
wish=read.csv("prediction.list",sep="\t")


#hy01 specific
####use these modified files with the variable name column cahnged so it matches the other files 
virhostmatcher=read.csv("hy01_virhostmatcher_gordonia_paper.melted.upset.csv")
phirbo=read.csv("hy01_phirbo_gordonia_paper.list.matrix.melted.upset.csv")
wish=read.csv("hy01_wish_prediction.upset.list",sep="\t")

#allviruses together specific
####use these modified files with the variable name column cahnged so it matches the other files 
virhostmatcher=read.csv("virhostmatcher_allviruses.csv")
phirbo=read.csv("phirbo_allviruses.csv")
wish=read.csv("wish_allviruses.csv")
php=read.csv("php_allviruses.csv")


setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/WISH/formanuscript/SFP10/")
all=read.csv("SFP10_confirmatory.csv")

library(data.table)
virhostmatcher=all[all$tool == "VHM",]
phirbo=all[all$tool == "Phirbo",]
php=all[all$tool == "PHP",]
wish=all[all$tool == "WIsH",]

#only GMAS/GRUS
virhostmatcher=all[all$tool == "VHM" & all$virus %like% "G",]
phirbo=all[all$tool == "Phirbo" & all$virus %like% "G",]
php=all[all$tool == "PHP" & all$virus %like% "G",]
wish=all[all$tool == "WIsH" & all$virus %like% "G",]

#only ecoli virus
virhostmatcher=all[all$tool == "VHM" & all$virus %like% "KFS|SFP|HY",]
phirbo=all[all$tool == "Phirbo" & all$virus %like% "KFS|SFP|HY",]
php=all[all$tool == "PHP" & all$virus %like% "KFS|SFP|HY",]
wish=all[all$tool == "WIsH" & all$virus %like% "KFS|SFP|HY",]


library(dplyr)
virhostmatcher <- virhostmatcher %>%
  mutate(category = case_when(score < 0.175 & shape == 16 ~ 'TP',
                              score < 0.175 & shape == 15 ~ 'FP',
                              score > 0.175 & shape == 15 ~ 'TN',
                              score > 0.175 & shape == 16 ~ 'FN'))

phirbo <- phirbo %>%
  mutate(category = case_when(score > 0.2 & shape == 16 ~ 'TP',
                              score > 0.2 & shape == 15 ~ 'FP',
                              score < 0.2 & shape == 15 ~ 'TN',
                              score < 0.2 & shape == 16 ~ 'FN'))

wish <- wish %>%
  mutate(category = case_when(pvalue < 0.06 & shape == 16 ~ 'TP',
                              pvalue < 0.06 & shape == 15 ~ 'FP',
                              pvalue > 0.06 & shape == 15 ~ 'TN',
                              pvalue > 0.06 & shape == 16 ~ 'FN'))

php <- php %>%
  mutate(category = case_when(score > 1442 & shape == 16 ~ 'TP',
                              score > 1442 & shape == 15 ~ 'FP',
                              score < 1442 & shape == 15 ~ 'TN',
                              score < 1442 & shape == 16 ~ 'FN'))


#calculate sensitivity and specificity for each tool 
table(virhostmatcher$category)
table(phirbo$category)
table(wish$category)
table(php$category)



library(stringr)
library(ComplexHeatmap)
set1 <- str_c(virhostmatcher$virus,"_", virhostmatcher$host,"_type_",virhostmatcher$category)
set2 <- str_c(phirbo$virus,"_", phirbo$host,"_type_",phirbo$category)
set3 <- str_c(wish$virus,"_", wish$host,"_type_",wish$category)
set4 <- str_c(php$virus,"_", php$host,"_type_",php$category)


#lt=list(set2,set4,set1,set3)
lt=list(set3,set1,set2,set4)
set_matrix=list_to_matrix(lt)
qwerty=as.data.frame(set_matrix)
library(tibble)
qwerty <- tibble::rownames_to_column(qwerty, "VALUE")
library(dplyr)
library(tidyr)
qwerty = qwerty %>%
  separate(VALUE, c("foo", "bar"), "_type_")

colnames(qwerty) = c("virus-host interaction","Accuracy","WIsH","VirHostMatcher","PHP","Phirbo")
interactions = colnames(qwerty)[3:6]
library(ComplexUpset)
library(ggplot2)
ComplexUpset::upset(
  qwerty,
  interactions,
  name='Confirmatory Tools',
  sort_intersections_by=c('degree', 'cardinality'),
set_sizes=FALSE,
  intersections='all',
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=Accuracy)
    )
  ),
  width_ratio=0.1 
)


#####
###########################################
############################
######
library(stringr)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)
        
#exploratroy tools upset
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/hostrange/")
all=read.csv("exploratoryresults_upset.csv")

#have to have them seperated by tool then glue back together
cherry=all[all$tool == "CHERRY",]
vhmn=all[all$tool == "VHMN",]
set1 <- str_c(cherry$virus,"_", cherry$host,"_type_",cherry$category)
set2 <- str_c(vhmn$virus,"_", vhmn$host,"_type_",vhmn$category)
lt=list(set1,set2)
set_matrix=list_to_matrix(lt)
qwerty=as.data.frame(set_matrix)
qwerty <- tibble::rownames_to_column(qwerty, "VALUE")
qwerty = qwerty %>%
  separate(VALUE, c("foo", "bar"), "_type_")


colnames(qwerty) = c("virus-host interaction","Accuracy","CHERRY","VHMN")
interactions = colnames(qwerty)[3:4]
library(ComplexUpset)
library(ggplot2)
ComplexUpset::upset(
  qwerty,
  interactions,
  name='Exploratory Tools',
  sort_intersections_by=c('degree', 'cardinality'),
  set_sizes=FALSE,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=Accuracy)) + coord_cartesian(ylim=c(0, 25))
    ),
  width_ratio=0.1 
)

table(cherry$category)
table(vhmn$category)
table(hostg$category)
table(rafah$category)
table(vpfclass$category)
table(vhulk$category)


hostg=all[all$tool == "HostG",]
rafah=all[all$tool == "RaFAH",]
vhulk=all[all$tool == "vHULK",]
vpfclass=all[all$tool == "vpf-class",]

set3 <- str_c(hostg$virus,"_", hostg$host,"_type_",hostg$category)
set4 <- str_c(rafah$virus,"_", rafah$host,"_type_",rafah$category)
set5 <- str_c(vhulk$virus,"_", vhulk$host,"_type_",vhulk$category)
set6 <- str_c(vpfclass$virus,"_", vpfclass$host,"_type_",vpfclass$category)

lt=list(set3,set4,set5,set6)

lt=list(set6,set5,set3,set3)

set_matrix=list_to_matrix(lt)
qwerty=as.data.frame(set_matrix)
qwerty <- tibble::rownames_to_column(qwerty, "VALUE")
qwerty = qwerty %>%
  separate(VALUE, c("foo", "bar"), "_type_")

colnames(qwerty) = c("virus-host interaction","Accuracy","HostG","RaFAH","vHULK","vpf-class")
colnames(qwerty) = c("virus-host interaction","Accuracy","HostG","RaFAH","vHULK","vpf-class")

interactions = colnames(qwerty)[3:6]
library(ComplexUpset)
library(ggplot2)
ComplexUpset::upset(
  qwerty,
  interactions,
  name='Exploratory Tools',
  sort_intersections_by=c('degree', 'cardinality'),
  set_sizes=FALSE,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=Accuracy)) + coord_cartesian(ylim=c(0, 25))
  ),
  width_ratio=0.1 
)



