#This should be run using Rscript and takes a command-line
#argument of a directory to process *.featurecount files
#from... the directory being a sample name
#This gets its "args" from the script run-R-sum.sh
#in the form of:
#Rscript --verbose $WORKERS_DIR/sumCounts.R $PWD $SAMPLE [Pfam|Kegg] $UPROC_FEATSUM
#Where PWD is the input dir
#And UPROC_FEATSUM is the output dir
#Sample and Pfam|Kegg are identifiers to prepend to output file

library(plyr)

args <- commandArgs(trailingOnly = TRUE)

setwd(paste(args[1]))

files=dir(pattern="*.featurecount")

holder <- read.delim(files[1],stringsAsFactors = FALSE,header=F)
holder$V2 <- 0

for (i in 1:length(files)) {
  tmp <- read.delim(files[i],stringsAsFactors = FALSE,header=F)
  tmp1=merge(holder,tmp,by = "V1",all=T)
  tmp1[is.na(tmp1)]<-0
  ddply(tmp1,.(V1),summarize,V2=sum(V2.x,V2.y))->holder
}

setwd(paste(args[4]))
write.table(holder,paste(args[3],args[2],"total_counts.tab",sep="_"),row.names=F,quote = F,sep = "\t")
