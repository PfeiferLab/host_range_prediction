library(intervals)
source("asynt.R")

###TO RUN
#SELCT A BLAST = ASSIGN EACH RESULTS A GAP1 GAP2 OR GAP3
#CONDENSE THOSE GAP DFS INTO GAPLAP1 (1V2) AND GAPLAP2 (1V2 V 3)
#FRM THOSE GAPS GET THE CORRESPONDINH 10536 PROKKA OVERLAPS

#####
#X vs Y means X is acting as the reference genome
blast=read.csv("blastn_15036_v_15144.blastn", sep="\t",header = FALSE)
blast=read.csv("blastn_15036_v_2192.blastn", sep="\t",header = FALSE)
blast=read.csv("blastn_15036_v_2196.blastn", sep="\t",header = FALSE)

#blast=read.csv("blastn_15144_v_10536.blastn", sep="\t",header = FALSE)

colnames(blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast$sstart2=ifelse(blast$sstart > blast$send, blast$send, blast$sstart)
blast$send2=ifelse(blast$send > blast$sstart, blast$send, blast$sstart)


#and make the plot
#par(mar = c(2,0,2,0))
#plot.alignments(alignments)
#alignments <- subset(alignments, Rlen >= 100 & query=="contig30.1" & reference == "mxdp_29")
#alignments <- subset(alignments, MQ >= 60)


#write output gap files cross-reffed between 1 infectable coli and 3 uninfectables
#check for gaps
library(SVbyEye)
target.ranges <- GenomicRanges::GRanges(
  seqnames = "target",
  IRanges::IRanges(
    #    start = alignments$Rstart,
    #    end = alignments$Rend
    
    start = blast$qstart,
    end = blast$qend
  )
)

query.ranges <- GenomicRanges::GRanges(
  seqnames = "query",
  IRanges::IRanges(
    #    start = alignments$Qstart,
    #   end = alignments$end
    start = blast$sstart2,
    end = blast$send2
    
  )
)


#do this for each of the 10536's vs each
gaps=reportGaps(t.ranges = target.ranges, q.ranges = query.ranges)
gaps=as.data.frame(gaps)

#have to flip the target and query to get gaps relative to one or the other
gaps2=reportGaps(t.ranges = target.ranges, q.ranges = query.ranges)
gaps2=as.data.frame(gaps2)

#have to flip the target and query to get gaps relative to one or the other
gaps3=reportGaps(t.ranges = target.ranges, q.ranges = query.ranges)
gaps3=as.data.frame(gaps3)

##compare gaps 1 v 2 
library(IRanges)
query <- IRanges(
  start = gaps$start,
  end = gaps$end
)

subject <- IRanges(
  start = gaps2$start,
  end = gaps2$end
)

#compare gaps 1v2 v 3
gaplap1=as.data.frame(IRanges::overlapsRanges(query,subject))

library(IRanges)
query <- IRanges(
  start = gaplap1$start,
  end = gaplap1$end
)

subject <- IRanges(
  start = gaps3$start,
  end = gaps3$end
)

gaplap2=as.data.frame(IRanges::overlapsRanges(query,subject))


#crossreference gaps with prokka gff file
library(IRanges)
  prokka=read.csv("PROKKA_05242023_10536.gff", sep="\t", header=FALSE)
prokka1=prokka[prokka$V1 == "f33647f64e1f4181_1",]

query <- IRanges(
  start = prokka1$V4,
  end = prokka1$V5
)

subject <- IRanges(
  start = gaplap2$start,
  end = gaplap2$end
)

prokka1$overlap=IRanges::overlapsAny(query, subject)

#prokka1$overlap=IRanges::overlapsAny(subject, query)

#prokka2 entries are fewer than crossref entries
#may be multiple gaps in same prokka gene entry
#so therefore cant easily do "prokka2$width=crossref$width" to filter short overlaps
prokka2=prokka1[prokka1$overlap == "TRUE",]
#crossref function not suoer useful at this omoment
#crossref=as.data.frame(IRanges::overlapsRanges(subject, query))
prokka2=within(prokka2, V9<-data.frame(do.call('rbind', strsplit(as.character(V9), 'Name=', fixed=TRUE))))
prokka2=within(prokka2, V9$X2<-data.frame(do.call('rbind', strsplit(as.character(V9$X2), ';', fixed=TRUE))))

#list of actual gene names
prokka3=prokka2[!grepl("ID=", prokka2$V9$X2$X1),]
write.csv(prokka3, "143nonhomolgousout.csv")
