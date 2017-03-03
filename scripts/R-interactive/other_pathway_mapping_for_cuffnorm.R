setwd("/Users/Scott/cuffnorm-out")

library(tidyr)
library(xlsx)
library(plyr)

#setup####
with_pathways<-read.table("with_pathways.tab",header=T)
attach(with_pathways)

#PATHWAYS THAT ARE UP OVERALL####

#oxidative phosphorylation####

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
#write.xlsx(oxphos_ec_sums,"oxphos_ec_sums_with_aliases.xlsx",row.names = F)

#heatmap for oxphos ####

oxphos2 <- oxphos_ec_sums[,c("ec_number","combined_effect","Hhep_effect","smad_effect")]

row.names(oxphos2)<-oxphos2$ec_number

oxphos2$gene_name = c("ppk","nox1, ahpF","nuo","frd","flii/yscn","sdhA","atpA","sdhB","nuoF","ppa")

oxphos2 <- oxphos2[order(oxphos2$combined_effect),]

x <- data.matrix(oxphos2[,2:4])

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none",margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))

row.names(oxphos2)<-oxphos2$gene_name

x <- data.matrix(oxphos2[,2:4])

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none",margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))

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

# alanine, aspartate and glutamate metabolism ####

alanine<-with_pathways[pathway=="00250|Alanine, aspartate and glutamate metabolism",]
no_blanks<-alanine[alanine$product_name!="",]
alanine<-no_blanks
rm(no_blanks)

alanine_ec_sums<-rowsum(alanine[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=alanine$ec_number)
truthy<-apply(alanine_ec_sums,1,function(row) all(row != 0))
no_zeros<-alanine_ec_sums[truthy,]
alanine_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(alanine_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
alanine_ec_sums$ec_number=row.names(alanine_ec_sums)

#decorate ec_sums with gene and product names
alanine_ec_sums <- merge(alanine_ec_sums,alanine[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-alanine_ec_sums[!duplicated(alanine_ec_sums),]
alanine_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(alanine_ec_sums,aliases,c(product_name,gene),sep = ",")
alanine_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(alanine_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
alanine_ec_sums<-collapsed
rm(collapsed)

alanine_ec_sums$combined_effect<-log(alanine_ec_sums$`S-H+ (Combined)`/alanine_ec_sums$`S+H- (Control)`)
alanine_ec_sums$Hhep_effect<-log(alanine_ec_sums$`S+H+ (H. hepaticus only)`/alanine_ec_sums$`S+H- (Control)`)
alanine_ec_sums$smad_effect<-log(alanine_ec_sums$`S-H- (SMAD3 Knockout)`/alanine_ec_sums$`S+H- (Control)`)
alanine_ec_sums<-alanine_ec_sums[order(alanine_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(alanine_ec_sums,"alanine_ec_sums_with_aliases.xlsx",row.names = F)

# 00630|Glyoxylate and dicarboxylate metabolism ####

glyoxylate<-with_pathways[pathway=="00630|Glyoxylate and dicarboxylate metabolism",]
no_blanks<-glyoxylate[glyoxylate$product_name!="",]
glyoxylate<-no_blanks
rm(no_blanks)

glyoxylate_ec_sums<-rowsum(glyoxylate[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=glyoxylate$ec_number)
truthy<-apply(glyoxylate_ec_sums,1,function(row) all(row != 0))
no_zeros<-glyoxylate_ec_sums[truthy,]
glyoxylate_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(glyoxylate_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
glyoxylate_ec_sums$ec_number=row.names(glyoxylate_ec_sums)

#decorate ec_sums with gene and product names
glyoxylate_ec_sums <- merge(glyoxylate_ec_sums,glyoxylate[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-glyoxylate_ec_sums[!duplicated(glyoxylate_ec_sums),]
glyoxylate_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(glyoxylate_ec_sums,aliases,c(product_name,gene),sep = ",")
glyoxylate_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(glyoxylate_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
glyoxylate_ec_sums<-collapsed
rm(collapsed)

glyoxylate_ec_sums$combined_effect<-log(glyoxylate_ec_sums$`S-H+ (Combined)`/glyoxylate_ec_sums$`S+H- (Control)`)
glyoxylate_ec_sums$Hhep_effect<-log(glyoxylate_ec_sums$`S+H+ (H. hepaticus only)`/glyoxylate_ec_sums$`S+H- (Control)`)
glyoxylate_ec_sums$smad_effect<-log(glyoxylate_ec_sums$`S-H- (SMAD3 Knockout)`/glyoxylate_ec_sums$`S+H- (Control)`)
glyoxylate_ec_sums<-glyoxylate_ec_sums[order(glyoxylate_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(glyoxylate_ec_sums,"glyoxylate_ec_sums_with_aliases.xlsx",row.names = F)

# 00680|Methane metabolism ####

methane<-with_pathways[pathway=="00680|Methane metabolism",]
no_blanks<-methane[methane$product_name!="",]
methane<-no_blanks
rm(no_blanks)

methane_ec_sums<-rowsum(methane[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=methane$ec_number)
truthy<-apply(methane_ec_sums,1,function(row) all(row != 0))
no_zeros<-methane_ec_sums[truthy,]
methane_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(methane_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
methane_ec_sums$ec_number=row.names(methane_ec_sums)

#decorate ec_sums with gene and product names
methane_ec_sums <- merge(methane_ec_sums,methane[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-methane_ec_sums[!duplicated(methane_ec_sums),]
methane_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(methane_ec_sums,aliases,c(product_name,gene),sep = ",")
methane_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(methane_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
methane_ec_sums<-collapsed
rm(collapsed)

methane_ec_sums$combined_effect<-log(methane_ec_sums$`S-H+ (Combined)`/methane_ec_sums$`S+H- (Control)`)
methane_ec_sums$Hhep_effect<-log(methane_ec_sums$`S+H+ (H. hepaticus only)`/methane_ec_sums$`S+H- (Control)`)
methane_ec_sums$smad_effect<-log(methane_ec_sums$`S-H- (SMAD3 Knockout)`/methane_ec_sums$`S+H- (Control)`)
methane_ec_sums<-methane_ec_sums[order(methane_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(methane_ec_sums,"methane_ec_sums_with_aliases.xlsx",row.names = F)


#PATHWAYS THAT ARE DOWN OVERALL####

#Peptidoglycan biosynthesis####

peptid<-with_pathways[pathway=="00550|Peptidoglycan biosynthesis",]
no_blanks<-peptid[peptid$product_name!="",]
peptid<-no_blanks
rm(no_blanks)

peptid_ec_sums<-rowsum(peptid[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=peptid$ec_number)
truthy<-apply(peptid_ec_sums,1,function(row) all(row != 0))
no_zeros<-peptid_ec_sums[truthy,]
peptid_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(peptid_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
peptid_ec_sums$ec_number=row.names(peptid_ec_sums)

#decorate ec_sums with gene and product names
peptid_ec_sums <- merge(peptid_ec_sums,peptid[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-peptid_ec_sums[!duplicated(peptid_ec_sums),]
peptid_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(peptid_ec_sums,aliases,c(product_name,gene),sep = ",")
peptid_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(peptid_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
peptid_ec_sums<-collapsed
rm(collapsed)

peptid_ec_sums$combined_effect<-log(peptid_ec_sums$`S-H+ (Combined)`/peptid_ec_sums$`S+H- (Control)`)
peptid_ec_sums$Hhep_effect<-log(peptid_ec_sums$`S+H+ (H. hepaticus only)`/peptid_ec_sums$`S+H- (Control)`)
peptid_ec_sums$smad_effect<-log(peptid_ec_sums$`S-H- (SMAD3 Knockout)`/peptid_ec_sums$`S+H- (Control)`)
peptid_ec_sums<-peptid_ec_sums[order(peptid_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(peptid_ec_sums,"peptid_ec_sums_with_aliases.xlsx",row.names = F)

#00290|Valine, leucine and isoleucine biosynthesis####

valine<-with_pathways[pathway=="00290|Valine, leucine and isoleucine biosynthesis",]
no_blanks<-valine[valine$product_name!="",]
valine<-no_blanks
rm(no_blanks)

valine_ec_sums<-rowsum(valine[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=valine$ec_number)
truthy<-apply(valine_ec_sums,1,function(row) all(row != 0))
no_zeros<-valine_ec_sums[truthy,]
valine_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(valine_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
valine_ec_sums$ec_number=row.names(valine_ec_sums)

#decorate ec_sums with gene and product names
valine_ec_sums <- merge(valine_ec_sums,valine[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-valine_ec_sums[!duplicated(valine_ec_sums),]
valine_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(valine_ec_sums,aliases,c(product_name,gene),sep = ",")
valine_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(valine_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
valine_ec_sums<-collapsed
rm(collapsed)

valine_ec_sums$combined_effect<-log(valine_ec_sums$`S-H+ (Combined)`/valine_ec_sums$`S+H- (Control)`)
valine_ec_sums$Hhep_effect<-log(valine_ec_sums$`S+H+ (H. hepaticus only)`/valine_ec_sums$`S+H- (Control)`)
valine_ec_sums$smad_effect<-log(valine_ec_sums$`S-H- (SMAD3 Knockout)`/valine_ec_sums$`S+H- (Control)`)
valine_ec_sums<-valine_ec_sums[order(valine_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(valine_ec_sums,"valine_ec_sums_with_aliases.xlsx",row.names = F)

#Fatty acid metabolism####

fatty<-with_pathways[pathway=="00071|Fatty acid metabolism",]
no_blanks<-fatty[fatty$product_name!="",]
fatty<-no_blanks
rm(no_blanks)

fatty_ec_sums<-rowsum(fatty[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=fatty$ec_number)
truthy<-apply(fatty_ec_sums,1,function(row) all(row != 0))
no_zeros<-fatty_ec_sums[truthy,]
fatty_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(fatty_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
fatty_ec_sums$ec_number=row.names(fatty_ec_sums)

#decorate ec_sums with gene and product names
fatty_ec_sums <- merge(fatty_ec_sums,fatty[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-fatty_ec_sums[!duplicated(fatty_ec_sums),]
fatty_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(fatty_ec_sums,aliases,c(product_name,gene),sep = ",")
fatty_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(fatty_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
fatty_ec_sums<-collapsed
rm(collapsed)

fatty_ec_sums$combined_effect<-log(fatty_ec_sums$`S-H+ (Combined)`/fatty_ec_sums$`S+H- (Control)`)
fatty_ec_sums$Hhep_effect<-log(fatty_ec_sums$`S+H+ (H. hepaticus only)`/fatty_ec_sums$`S+H- (Control)`)
fatty_ec_sums$smad_effect<-log(fatty_ec_sums$`S-H- (SMAD3 Knockout)`/fatty_ec_sums$`S+H- (Control)`)
fatty_ec_sums<-fatty_ec_sums[order(fatty_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(fatty_ec_sums,"fatty_ec_sums_with_aliases.xlsx",row.names = F)

#00260|Glycine, serine and threonine metabolism####

glycine<-with_pathways[pathway=="00260|Glycine, serine and threonine metabolism",]
no_blanks<-glycine[glycine$product_name!="",]
glycine<-no_blanks
rm(no_blanks)

glycine_ec_sums<-rowsum(glycine[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=glycine$ec_number)
truthy<-apply(glycine_ec_sums,1,function(row) all(row != 0))
no_zeros<-glycine_ec_sums[truthy,]
glycine_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(glycine_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
glycine_ec_sums$ec_number=row.names(glycine_ec_sums)

#decorate ec_sums with gene and product names
glycine_ec_sums <- merge(glycine_ec_sums,glycine[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-glycine_ec_sums[!duplicated(glycine_ec_sums),]
glycine_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(glycine_ec_sums,aliases,c(product_name,gene),sep = ",")
glycine_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(glycine_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
glycine_ec_sums<-collapsed
rm(collapsed)

glycine_ec_sums$combined_effect<-log(glycine_ec_sums$`S-H+ (Combined)`/glycine_ec_sums$`S+H- (Control)`)
glycine_ec_sums$Hhep_effect<-log(glycine_ec_sums$`S+H+ (H. hepaticus only)`/glycine_ec_sums$`S+H- (Control)`)
glycine_ec_sums$smad_effect<-log(glycine_ec_sums$`S-H- (SMAD3 Knockout)`/glycine_ec_sums$`S+H- (Control)`)
glycine_ec_sums<-glycine_ec_sums[order(glycine_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(glycine_ec_sums,"glycine_ec_sums_with_aliases.xlsx",row.names = F)

#00195|Photosynthesis####

photosynth<-with_pathways[pathway=="00195|Photosynthesis",]
no_blanks<-photosynth[photosynth$product_name!="",]
photosynth<-no_blanks
rm(no_blanks)

photosynth_ec_sums<-rowsum(photosynth[,c("S3_FPM","S4_FPM","S1_FPM","S2_FPM")],group=photosynth$ec_number)
truthy<-apply(photosynth_ec_sums,1,function(row) all(row != 0))
no_zeros<-photosynth_ec_sums[truthy,]
photosynth_ec_sums<-no_zeros
rm(truthy,no_zeros)

colnames(photosynth_ec_sums)=c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
photosynth_ec_sums$ec_number=row.names(photosynth_ec_sums)

#decorate ec_sums with gene and product names
photosynth_ec_sums <- merge(photosynth_ec_sums,photosynth[,c("product_name","gene","ec_number")],by="ec_number",all=F)
no_dups<-photosynth_ec_sums[!duplicated(photosynth_ec_sums),]
photosynth_ec_sums<-no_dups
rm(no_dups)

#combine gene and product_name columns
combined<-unite(photosynth_ec_sums,aliases,c(product_name,gene),sep = ",")
photosynth_ec_sums<-combined
rm(combined)

#grab together the multiple product names
collapsed<-ddply(photosynth_ec_sums, .(ec_number,`S+H- (Control)`,`S-H- (SMAD3 Knockout)`,`S+H+ (H. hepaticus only)`,`S-H+ (Combined)`), summarize, aliases=toString(aliases))
photosynth_ec_sums<-collapsed
rm(collapsed)

photosynth_ec_sums$combined_effect<-log(photosynth_ec_sums$`S-H+ (Combined)`/photosynth_ec_sums$`S+H- (Control)`)
photosynth_ec_sums$Hhep_effect<-log(photosynth_ec_sums$`S+H+ (H. hepaticus only)`/photosynth_ec_sums$`S+H- (Control)`)
photosynth_ec_sums$smad_effect<-log(photosynth_ec_sums$`S-H- (SMAD3 Knockout)`/photosynth_ec_sums$`S+H- (Control)`)
photosynth_ec_sums<-photosynth_ec_sums[order(photosynth_ec_sums$combined_effect , decreasing = T),]

#lets write some excel tables
write.xlsx(photosynth_ec_sums,"photosynth_ec_sums_with_aliases.xlsx",row.names = F)
