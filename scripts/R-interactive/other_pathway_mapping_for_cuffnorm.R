setwd("/Users/Scott/cuffnorm-out")

library(reshape2)
library(tidyr)
library(RColorBrewer)
library(xlsx)

#setup####
with_pathways<-read.table("with_pathways.tab",header=T)

#oxidative phosphorylation####

attach(with_pathways)

oxphos<-with_pathways[pathway=="00190|Oxidative phosphorylation",]
no_blanks<-oxphos[oxphos$product_name!="",]
oxphos<-no_blanks
rm(no_blanks)

oxphos_product_sums<-rowsum(oxphos[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=oxphos$product_name)
truthy<-apply(oxphos_product_sums,1,function(row) all(row != 0))
no_zeros<-oxphos_product_sums[truthy,]
oxphos_product_sums<-no_zeros
rm(truthy,no_zeros)

colnames(oxphos_product_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")

oxphos_product_sums$combined_effect<-log(oxphos_product_sums$`S-H+ (Combined)`/oxphos_product_sums$`S+H- (Control)`)
oxphos_product_sums$Hhep_effect<-log(oxphos_product_sums$`S+H+ (H. hepaticus only)`/oxphos_product_sums$`S+H- (Control)`)
oxphos_product_sums$smad_effect<-log(oxphos_product_sums$`S-H- (SMAD3 Knockout)`/oxphos_product_sums$`S+H- (Control)`)
oxphos_product_sums<-oxphos_product_sums[order(oxphos_product_sums$combined_effect , decreasing = T),]

#welp, the top changed product is atp synthase, both positive and negative! ha!

#what are the actual species for atpb
# genome_to_feature <- read.delim("genome-name_to_refseq-locus-tag")
#
# atpb<-with_pathways[product_name=="atp synthase beta chain" & pathway=="00190|Oxidative phosphorylation",]
# important <- merge(x=atpb,y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all = F)
# #and for the lowest f0f1 atp synthase subunit alpha
# atpa<-with_pathways[product_name=="atp synthase alpha chain" & pathway=="00190|Oxidative phosphorylation",]
# important2 <- merge(x=atpa,y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all = F)

#lets write some excel tables
# write.xlsx(oxphos_product_sums,"oxphos_product_sums.xlsx")
# write.xlsx(important,"atp synthase beta chain downregulated.xlsx")
# write.xlsx(important2,"atp synthatse alpha chain upregualted.xlsx")

# benzoate degredation ####

benzoate<-with_pathways[pathway=="00362|Benzoate degradation via hydroxylation",]
no_blanks<-benzoate[benzoate$product_name!="",]
benzoate<-no_blanks
rm(no_blanks)

benzoate_product_sums<-rowsum(benzoate[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=benzoate$product_name)
truthy<-apply(benzoate_product_sums,1,function(row) all(row != 0))
no_zeros<-benzoate_product_sums[truthy,]
benzoate_product_sums<-no_zeros
rm(truthy,no_zeros)

colnames(benzoate_product_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")

benzoate_product_sums$combined_effect<-log(benzoate_product_sums$`S-H+ (Combined)`/benzoate_product_sums$`S+H- (Control)`)
benzoate_product_sums$Hhep_effect<-log(benzoate_product_sums$`S+H+ (H. hepaticus only)`/benzoate_product_sums$`S+H- (Control)`)
benzoate_product_sums$smad_effect<-log(benzoate_product_sums$`S-H- (SMAD3 Knockout)`/benzoate_product_sums$`S+H- (Control)`)
benzoate_product_sums<-benzoate_product_sums[order(benzoate_product_sums$combined_effect , decreasing = T),]



#lets write some excel tables
write.xlsx(benzoate_product_sums,"benzoate_product_sums.xlsx")
