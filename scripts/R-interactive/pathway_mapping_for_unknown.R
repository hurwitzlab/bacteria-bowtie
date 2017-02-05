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

patric_annotation <- read.delim("../unknown-cuffnorm-out/product_to_pathway.tab")
patric_annotation$product <- tolower(patric_annotation$product)
with_pathways<-merge(x=sum_by_product_name,y=patric_annotation,by.x="product",by.y="product")
with_pathways_dups_removed<-with_pathways[!duplicated(with_pathways),]
with_pathways<-with_pathways_dups_removed
rm(with_pathways_dups_removed)

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
more_shortened<-shortened[1:30,]

setwd("~/bacteria-bowtie/scripts/R-interactive/")
write.table(more_shortened,"unknown_sum_by_kegg_pathway.tab", sep = "\t", quote = T,row.names = F)

system("source ~/.bash_profile && ./bubble.sh unknown_sum_by_kegg_pathway.tab unknownBubble")

system("cp unknownBubble.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")

more_shortened$unknown_effect<-log(more_shortened$`S-H+ (unknown)`/more_shortened$`S+H- (Control)`)
more_shortened$Hhep_effect<-log(more_shortened$`S+H+ (H. hepaticus only)`/more_shortened$`S+H- (Control)`)
more_shortened$smad_effect<-log(more_shortened$`S-H- (SMAD3 Knockout)`/more_shortened$`S+H- (Control)`)
more_shortened<-more_shortened[order(more_shortened$unknown_effect , decreasing = T),]

#tom's request
more_shortened$control_effect<-0

effects <- data.matrix(more_shortened[,6:9])

row.names(effects)<-more_shortened$Name

x=effects

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

#a heatmap, cuz why not!
heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("combined", "H. hepaticus","SMAD3-KO","Nothing"))

new_bubble_source <- more_shortened[,1:5]

write.table(new_bubble_source,"sum_by_kegg_pathway_ordered_by_combined_effect.tab", sep = "\t", quote = T,row.names = F)

system("source ~/.bash_profile && ./bubble.sh sum_by_kegg_pathway_ordered_by_combined_effect.tab unknownBubble_orderedbyeffect")

system("cp unknownBubble_orderedbyeffect.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")
