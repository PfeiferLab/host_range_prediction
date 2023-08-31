library(intervals)
source("asynt.R")

###TO RUN
#SELCT A BLAST = ASSIGN EACH RESULTS A prokka3 prokka4 or prokka5 
#find the shared gene names from those prokkas

#####
#X vs Y means X is acting as the reference genome 
blast=read.csv("blastn_2192_v_10536.blastn", sep="\t",header = FALSE)
blast=read.csv("blastn_2196_v_10536.blastn", sep="\t",header = FALSE)
blast=read.csv("blastn_15144_v_10536.blastn", sep="\t",header = FALSE)

colnames(blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast$sstart2=ifelse(blast$sstart > blast$send, blast$send, blast$sstart)
blast$send2=ifelse(blast$send > blast$sstart, blast$send, blast$sstart)

#check for gaps
library(SVbyEye)
target.ranges <- GenomicRanges::GRanges(
  seqnames = "target",
  IRanges::IRanges(
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


#gaps f33 are clear in the circo plot but not gaps from b70 
#f33 is reference 10536 which is the one it CAN infect 

#have to flip the target and query to get gaps relative to one or the other
gaps=reportGaps(t.ranges = target.ranges, q.ranges = query.ranges)
gaps=as.data.frame(gaps)


#crossreference gaps with prokka gff file
library(IRanges)
prokka=read.csv("PROKKA_05242023_2192.gff", sep="\t", header=FALSE)
prokka1=prokka[prokka$V1 == "4d52a7c47b42404f_1",]

prokka=read.csv("PROKKA_05242023_2196.gff", sep="\t", header=FALSE)
prokka1=prokka[prokka$V1 == "0795399960f546b3_1",]

prokka=read.csv("PROKKA_05242023_15144.gff", sep="\t", header=FALSE)
prokka1=prokka[prokka$V1 == "b70a1a81a2744334_1",]

query <- IRanges(
  start = prokka1$V4,
  end = prokka1$V5
)

subject <- IRanges(
  start = gaps$start,
  end = gaps$end
)

prokka1$overlap=IRanges::overlapsAny(query, subject)

#prokka1$overlap=IRanges::overlapsAny(subject, query)

#prokka2 entries are fewer than crossref entries
#may be multiple gaps in same prokka gene entry 
#so therefore cant easily do "prokka2$width=crossref$width" to filter short overlaps
prokka2=prokka1[prokka1$overlap == "TRUE",]
#crossref function not useful at this omoment 
#crossref=as.data.frame(IRanges::overlapsRanges(subject, query))
prokka2=within(prokka2, V9<-data.frame(do.call('rbind', strsplit(as.character(V9), 'Name=', fixed=TRUE))))
prokka2=within(prokka2, V9$X2<-data.frame(do.call('rbind', strsplit(as.character(V9$X2), ';', fixed=TRUE))))

#list of actual gene names 
prokka3=prokka2[!grepl("ID=", prokka2$V9$X2$X1),]
#get for each X v. 10536 the genes and see what gene name overlap there is betweeen the 3 
prokka4=prokka2[!grepl("ID=", prokka2$V9$X2$X1),]
prokka5=prokka2[!grepl("ID=", prokka2$V9$X2$X1),]

combine1=intersect(prokka3$V9$X2$X1, prokka4$V9$X2$X1) 
combine2=intersect(combine1, prokka5$V9$X2$X1)


#pull out what the functions are from any prokka file doesnt amtter
matches=prokka3[grep(paste(combine2, collapse='|'), prokka3$V9$X2$X1, ignore.case=TRUE),]

write.csv(matches, "supplementarytableX.csv")
