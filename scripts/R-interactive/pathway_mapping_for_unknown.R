#DONE
#these commands make id_to_gene, id_to_product and id_to_ec_number
#setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")
#system(command = "./cuffnorm-pathways.sh")
#DONE

library(reshape2)
library(tidyr)

setwd("/Users/Scott/unknown-cuffnorm-out/")

#setup####
filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
filtered_annotated <- filtered_annotated[,2:8]
attach(filtered_annotated)

sum_by_product_name <- read.csv("sum_by_product_name.csv")
colnames(sum_by_product_name)[1]<-"product"
attach(sum_by_product_name)

sum_by_product_name<-sum_by_product_name[,1:5]
sum_by_product_name$product <- tolower(sum_by_product_name$product)

#now we will attempt to add pathways####

patric_annotation <- read.delim("../combined-cuffnorm-out/product_to_pathway.tab")
patric_annotation$product <- tolower(patric_annotation$product)
with_pathways<-merge(x=sum_by_product_name,y=patric_annotation,by.x="product",by.y="product")

####Start here again####
#dont want blanks
with_pathways<-with_pathways[grep(".+",with_pathways$pathway),]
#separate into individual pathways
with_pathways<-separate_rows(with_pathways, pathway, sep = ";")
sum_by_kegg_pathway<-rowsum(with_pathways[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways$pathway)
sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])

#formatting for stupid bubble plot, which i bet could be done in R
sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)

shortened<-sum_by_kegg_pathway[sum_by_kegg_pathway$sum>mean(sum_by_kegg_pathway$sum),]
shortened<-shortened[order(shortened$sum , decreasing = T),]
shortened<-shortened[,c("Name","S3_FPM","S4_FPM","S1_FPM","S2_FPM")]
colnames(shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (unknown)")

setwd("~/bacteria-bowtie/scripts/R-interactive/")
write.table(shortened,"unknown_sum_by_kegg_pathway_above_mean.tab", sep = "\t", quote = T,row.names = F)

shortened<-read.table("unknown_sum_by_kegg_pathway_above_mean.tab",header = T)
more_shortened<-shortened[1:30,]
colnames(more_shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (unknown)")

setwd("~/bacteria-bowtie/scripts/R-interactive/")
write.table(more_shortened,"unknown_sum_by_kegg_pathway_above_mean.tab", sep = "\t", quote = T,row.names = F)

#Have to run this in an external shell cuz ... perl... grumble
#system("./bubble.sh unknown_sum_by_kegg_pathway_above_mean.tab unknownBubble")

system("cp unknownBubble.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")

more_shortened$unknown_effect<-log(more_shortened$`S-H+ (unknown)`/more_shortened$`S+H- (Control)`)
more_shortened<-more_shortened[order(more_shortened$unknown_effect , decreasing = T),]

#top_sixty_five<-sum_by_kegg_pathway[order(sum_by_kegg_pathway$sum,decreasing = T)[1:65],]
#top_sixty_five<-top_sixty_five[,c("Name","S1_FPM","S2_FPM","S3_FPM","S4_FPM")]

#write.table(top_sixty_five,"~/tophat-bacteria/scripts/R-interactive/sum_by_kegg_pathway_top_sixty_five.tab", sep = "\t", quote = T,row.names = F)

