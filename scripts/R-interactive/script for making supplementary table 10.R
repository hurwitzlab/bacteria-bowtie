setwd("/Users/Scott/cuffnorm-out")

library(tidyr)
library(xlsx)
library(plyr)

#setup####
with_pathways<-read.table("with_pathways.tab",header=T)
attach(with_pathways)

sum_by_kegg_pathway<-rowsum(with_pathways[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways$pathway)
sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])
sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)

everyting<-sum_by_kegg_pathway
everyting<-everyting[order(everyting$sum , decreasing = T),]

everyting<-everyting[,c("Name","S3_FPM","S4_FPM","S1_FPM","S2_FPM")]
colnames(everyting)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")

everyting$combined_effect<-log(everyting$`S-H+ (Combined)`/everyting$`S+H- (Control)`)
everyting$Hhep_effect<-log(everyting$`S+H+ (H. hepaticus only)`/everyting$`S+H- (Control)`)
everyting$smad_effect<-log(everyting$`S-H- (SMAD3 Knockout)`/everyting$`S+H- (Control)`)
everyting<-everyting[order(everyting$combined_effect , decreasing = T),]

write.xlsx(everyting,row.names = F,col.names = T,file = "140_pathways.xlsx")
