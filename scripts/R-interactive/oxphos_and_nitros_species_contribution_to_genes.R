setwd("/Users/Scott/cuffnorm-out")

library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggiraph)

#setup####
with_pathways<-read.table("with_pathways.tab",header=T)
attach(with_pathways)
genome_to_feature <- read.delim("genome-name_to_refseq-locus-tag",colClasses = "character")

#oxphos####
oxphos<-with_pathways[pathway=="00190|Oxidative phosphorylation",]
no_blanks<-oxphos[oxphos$product_name!="",]
oxphos<-no_blanks
rm(no_blanks)

detach(with_pathways)
attach(oxphos)

important<-oxphos[ec_number %in% c("2.7.4.1","1.6.99.3","1.6.5.3"),]
important2 <- merge(x=important[,c(1:5,9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag")

colnames(important2)[2:5]<-c("S+H+ (H. hepaticus)","S-H+ (Combined)","S+H- (Control)","S-H- (Smad3-/-)")

important2$genome_name<-as.character(important2$genome_name)
melted<-melt(important2[,c(2:7)])
order<-c("S+H- (Control)","S-H- (Smad3-/-)","S+H+ (H. hepaticus)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("ec_number","genome","Sample","RNA count")

#going down to just 1.6.5.3 and 2.7.4.1
melted <- melted[melted$ec_number %in% c("1.6.5.3","2.7.4.1"),]
melted$ec_number = gsub('1.6.5.3','[EC:1.6.5.3]',melted$ec_number)
melted$ec_number = gsub('2.7.4.1','[EC:2.7.4.1]',melted$ec_number)
melted <- melted[melted$`RNA count` != 0,]

melted <- with(melted, melted[order(ec_number, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

#interactive
plot.bar <- ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))




#nitros####
detach(oxphos)
attach(with_pathways)

nitros<-with_pathways[pathway=="00910|Nitrogen metabolism",]
no_blanks<-nitros[nitros$product_name!="",]
nitros<-no_blanks
rm(no_blanks)

detach(with_pathways)
attach(nitros)

#important<-nitros[ec_number %in% c("4.3.1.1","6.3.1.1"),]
important<-nitros[ec_number %in% c("1.4.1.13","1.6.5.3"),]
important2 <- merge(x=important[,c(1:5,9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag")

colnames(important2)[2:5]<-c("S+H+ (H. hepaticus)","S-H+ (Combined)","S+H- (Control)","S-H- (Smad3-/-)")

important2$genome_name<-as.character(important2$genome_name)
melted<-melt(important2[,c(2:7)])
order<-c("S+H- (Control)","S-H- (Smad3-/-)","S+H+ (H. hepaticus)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("ec_number","genome","Sample","count")

melted <- with(melted, melted[order(ec_number, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=count, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

#interactive
plot.bar <- ggplot(data=melted, aes(x=Sample, y=count, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~ec_number) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))

