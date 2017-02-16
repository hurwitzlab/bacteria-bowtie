setwd("/Users/Scott/cuffnorm-out")

library(tidyr)
library(xlsx)
library(plyr)

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
write.xlsx(oxphos_product_sums,"oxphos_product_sums.xlsx")
# write.xlsx(important,"atp synthase beta chain downregulated.xlsx")
# write.xlsx(important2,"atp synthatse alpha chain upregualted.xlsx")

#oxphos ec sums ####

oxphos_ec_sums<-rowsum(oxphos[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=oxphos$ec_number)
truthy<-apply(oxphos_ec_sums,1,function(row) all(row != 0))
no_zeros<-oxphos_ec_sums[truthy,]
oxphos_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(oxphos_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
oxphos_ec_sums$ec_number=row.names(oxphos_ec_sums)

#decorate ec_sums with gene and product names
oxphos_ec_sums <- merge(oxphos_ec_sums,oxphos[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-oxphos_ec_sums[!duplicated(oxphos_ec_sums),]
oxphos_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(oxphos_ec_sums,aliases,c(product_name,gene),sep = ",")
oxphos_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(oxphos_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
oxphos_ec_sums<-collapsed
rm(collapsed)

oxphos_ec_sums$combined_effect<-log(oxphos_ec_sums$`S-H+ (Combined)`/oxphos_ec_sums$`S+H- (Control)`)
oxphos_ec_sums$Hhep_effect<-log(oxphos_ec_sums$`S+H+ (H. hepaticus only)`/oxphos_ec_sums$`S+H- (Control)`)
oxphos_ec_sums$smad_effect<-log(oxphos_ec_sums$`S-H- (SMAD3 Knockout)`/oxphos_ec_sums$`S+H- (Control)`)
oxphos_ec_sums<-oxphos_ec_sums[order(oxphos_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(oxphos_ec_sums,"oxphos_ec_sums_with_aliases.xlsx",row.names = F)

#nitrogen metabolism ####

nitros<-with_pathways[pathway=="00910|Nitrogen metabolism",]
no_blanks<-nitros[nitros$product_name!="",]
nitros<-no_blanks
rm(no_blanks)

nitros_ec_sums<-rowsum(nitros[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=nitros$ec_number)
truthy<-apply(nitros_ec_sums,1,function(row) all(row != 0))
no_zeros<-nitros_ec_sums[truthy,]
nitros_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(nitros_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
nitros_ec_sums$ec_number=row.names(nitros_ec_sums)

#decorate ec_sums with gene and product names
nitros_ec_sums <- merge(nitros_ec_sums,nitros[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-nitros_ec_sums[!duplicated(nitros_ec_sums),]
nitros_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(nitros_ec_sums,aliases,c(product_name,gene),sep = ",")
nitros_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(nitros_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
nitros_ec_sums<-collapsed
rm(collapsed)

nitros_ec_sums$combined_effect<-log(nitros_ec_sums$`S-H+ (Combined)`/nitros_ec_sums$`S+H- (Control)`)
nitros_ec_sums$Hhep_effect<-log(nitros_ec_sums$`S+H+ (H. hepaticus only)`/nitros_ec_sums$`S+H- (Control)`)
nitros_ec_sums$smad_effect<-log(nitros_ec_sums$`S-H- (SMAD3 Knockout)`/nitros_ec_sums$`S+H- (Control)`)
nitros_ec_sums<-nitros_ec_sums[order(nitros_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(nitros_ec_sums,"nitros_ec_sums_with_aliases.xlsx",row.names = F)

#citrate tca cycle####

citrate<-with_pathways[pathway=="00020|Citrate cycle (TCA cycle)",]
no_blanks<-citrate[citrate$product_name!="",]
citrate<-no_blanks
rm(no_blanks)

citrate_ec_sums<-rowsum(citrate[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=citrate$ec_number)
truthy<-apply(citrate_ec_sums,1,function(row) all(row != 0))
no_zeros<-citrate_ec_sums[truthy,]
citrate_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(citrate_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
citrate_ec_sums$ec_number=row.names(citrate_ec_sums)

#decorate ec_sums with gene and product names
citrate_ec_sums <- merge(citrate_ec_sums,citrate[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-citrate_ec_sums[!duplicated(citrate_ec_sums),]
citrate_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(citrate_ec_sums,aliases,c(product_name,gene),sep = ",")
citrate_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(citrate_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
citrate_ec_sums<-collapsed
rm(collapsed)

citrate_ec_sums$combined_effect<-log(citrate_ec_sums$`S-H+ (Combined)`/citrate_ec_sums$`S+H- (Control)`)
citrate_ec_sums$Hhep_effect<-log(citrate_ec_sums$`S+H+ (H. hepaticus only)`/citrate_ec_sums$`S+H- (Control)`)
citrate_ec_sums$smad_effect<-log(citrate_ec_sums$`S-H- (SMAD3 Knockout)`/citrate_ec_sums$`S+H- (Control)`)
citrate_ec_sums<-citrate_ec_sums[order(citrate_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(citrate_ec_sums,"citrate_ec_sums_with_aliases.xlsx",row.names = F)


