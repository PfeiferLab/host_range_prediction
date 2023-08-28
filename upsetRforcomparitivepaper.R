## Written by: Abby Howell
## Modified by: Cyril Versoza
## Date: August 28, 2023

# Load necessary libraries:
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

# Assign variables:
virhostmatcher=read.csv("virhostmatcher_gordonia_paper.melted.csv")
phirbo=read.csv("phirbo_gordonia_paper.list.matrix.melted.csv")
wish=read.csv("prediction.list",sep="\t")


# HY01 specific
#### Use these modified files with the variable name column changed, so it matches the other files 
virhostmatcher=read.csv("hy01_virhostmatcher_gordonia_paper.melted.upset.csv")
phirbo=read.csv("hy01_phirbo_gordonia_paper.list.matrix.melted.upset.csv")
wish=read.csv("hy01_wish_prediction.upset.list",sep="\t")

# All viruses together; specific
#### Use these modified files with the variable name column changed, so it matches the other files 
virhostmatcher=read.csv("virhostmatcher_allviruses.csv")
phirbo=read.csv("phirbo_allviruses.csv")
wish=read.csv("wish_allviruses.csv")
php=read.csv("php_allviruses.csv")

all=read.csv("SFP10_confirmatory.csv")

virhostmatcher=all[all$tool == "VHM",]
phirbo=all[all$tool == "Phirbo",]
php=all[all$tool == "PHP",]
wish=all[all$tool == "WIsH",]

table(all$tool)

# Only GMAs/GRUs:
virhostmatcher=all[all$tool == "VHM" & all$virus %like% "G",]
phirbo=all[all$tool == "Phirbo" & all$virus %like% "G",]
php=all[all$tool == "PHP" & all$virus %like% "G",]
wish=all[all$tool == "WIsH" & all$virus %like% "G",]

# Only E coli virus:
virhostmatcher=all[all$tool == "VHM" & all$virus %like% "KFS|SFP|HY",]
phirbo=all[all$tool == "Phirbo" & all$virus %like% "KFS|SFP|HY",]
php=all[all$tool == "PHP" & all$virus %like% "KFS|SFP|HY",]
wish=all[all$tool == "WIsH" & all$virus %like% "KFS|SFP|HY",]

virhostmatcher$score=as.numeric(virhostmatcher$score)
virhostmatcher <- virhostmatcher %>%
  mutate(category = case_when(score > 0.175 & shape == 16 ~ 'TP',
                              score > 0.175 & shape == 15 ~ 'FP',
                              score < 0.175 & shape == 15 ~ 'TN',
                              score < 0.175 & shape == 16 ~ 'FN'))

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

phirbo=read.csv("SFP10_confirmatory_phirbo.csv")

# Calculate sensitivity and specificity for each tool 
table(virhostmatcher$category)
table(phirbo$category)
table(wish$category)
table(php$category)

set1 <- str_c(phirbo$virus,"_", phirbo$host,"_type_",phirbo$category)
set2 <- str_c(php$virus,"_", php$host,"_type_",php$category)
set3 <- str_c(virhostmatcher$virus,"_", virhostmatcher$host,"_type_",virhostmatcher$category)
set4 <- str_c(wish$virus,"_", wish$host,"_type_",wish$category)


lt=list(set4,set3,set2,set1)
set_matrix=list_to_matrix(lt)
qwerty=as.data.frame(set_matrix)
qwerty <- tibble::rownames_to_column(qwerty, "VALUE")
qwerty = qwerty %>%
  separate(VALUE, c("foo", "bar"), "_type_")


# Write threshold passing table:
passing=virhostmatcher[virhostmatcher$category == "TP" | virhostmatcher$category == "FP",]
passing2=php[php$category == "TP" | php$category == "FP",]
passing3=wish[wish$category == "TP" | wish$category == "FP",]
passing4=phirbo[phirbo$category == "TP" | phirbo$category == "FP",]


write.table( passing,
             file = "thresholdpassing.txt", 
             append = T,
             sep = " ")

write.table( passing2,
             file = "thresholdpassing.txt", 
             append = T,
             sep = " ")

write.table( passing3,
             file = "thresholdpassing.txt", 
             append = T,
             sep = " ")

write.table( passing4,
             file = "thresholdpassing.txt", 
             append = T,
             sep = " ")



colnames(qwerty) = c("virus-host interaction","Accuracy","WIsH","VirHostMatcher","PHP","Phirbo")
interactions = colnames(qwerty)[3:6]
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
  width_ratio=0.1,
sort_sets=FALSE

)

###########################


# Exploratroy tools -- upset
all=read.csv("exploratoryresults_upset.csv")

# Have to have them seperated by tool then glue back together
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
  width_ratio=0.1,
  sort_sets=FALSE 
)

hostg=all[all$tool == "HostG",]
rafah=all[all$tool == "RaFAH",]
vhulk=all[all$tool == "vHULK",]
vpfclass=all[all$tool == "vpf-class",]

table(cherry$category)
table(vhmn$category)
table(hostg$category)
table(rafah$category)
table(vpfclass$category)
table(vhulk$category)


set3 <- str_c(hostg$virus,"_", hostg$host,"_type_",hostg$category)
set4 <- str_c(rafah$virus,"_", rafah$host,"_type_",rafah$category)
set5 <- str_c(vhulk$virus,"_", vhulk$host,"_type_",vhulk$category)
set6 <- str_c(vpfclass$virus,"_", vpfclass$host,"_type_",vpfclass$category)

lt=list(set6,set5,set4,set3)

set_matrix=list_to_matrix(lt)
qwerty=as.data.frame(set_matrix)
qwerty <- tibble::rownames_to_column(qwerty, "VALUE")
qwerty = qwerty %>%
  separate(VALUE, c("foo", "bar"), "_type_")

colnames(qwerty) = c("virus-host interaction","Accuracy","vpf-class","vHULK","RaFAH","HostG")

interactions = colnames(qwerty)[3:6]
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
  width_ratio=0.1,
  sort_sets=FALSE
)

data=read.csv("vpfclass_all_host_genus.csv", sep="\t")
high=data[as.numeric(data$membership_ratio) >= 0.3 & data$confidence_score >= 0.5 ,]
