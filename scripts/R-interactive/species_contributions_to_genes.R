setwd("/Users/Scott/cuffnorm-out")

library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggiraph)

#probably don't need these yet
#library(KEGGgraph)
#library(KEGGREST)

#this might work later if I really wanted to work at it
#but too much work for now
#source("/Users/Scott/tophat-bacteria/scripts/R-interactive/uglyMerge.R")

#setup####
filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
filtered_annotated <- filtered_annotated[,2:9]
sum_by_product_name <- read.csv("sum_by_product_name.csv")
colnames(sum_by_product_name)[1]<-"product"
sum_by_product_name<-sum_by_product_name[,1:5]
sum_by_product_name$product <- tolower(sum_by_product_name$product)
sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
colnames(sum_by_gene_name)[1]<-"gene"
sum_by_gene_name<-sum_by_gene_name[,1:5]
sum_by_gene_name$gene <- tolower(sum_by_gene_name$gene)

#Got all.PATRIC.cds.tab from the server and cut down to just genome names and refseq locus tags
#system("cut -f2,7 all.PATRIC.cds.tab > genome-name_to_refseq-locus-tag")
genome_to_feature <- read.delim("genome-name_to_refseq-locus-tag")

#Let's try polyamines now####

polyamines<-read.table("polyamine-search-list",sep="\t",comment.char = "#")
polyamines<-data.frame(tolower(polyamines[,1]))
colnames(polyamines)<-"V1"
polyamines<-unique(polyamines)
poly_product<-data.frame(product=polyamines$V1[1:88])
poly_genes<-data.frame(gene=polyamines$V1[89:203])

just_poly_products<-merge(x=sum_by_product_name,y=poly_product,all=F)
just_poly_genes<-merge(x=sum_by_gene_name,y=poly_genes,all=F)

#getting more annotation for poly####
poly_annot<-read.table("polyamine_list_annotation",header = T,sep = ";",strip.white = T)
poly_annot$product <- tolower(poly_annot$product)
just_poly_products_annot <- merge(x=poly_annot,y=just_poly_products)
attach(just_poly_products_annot)

important<-just_poly_products_annot[ecnumber %in% c("[EC:3.5.1.53]","[EC:4.1.1.96]"),]
#doing those three that are in a row before putrescine
#important<-just_poly_products_annot[ecnumber %in% c("[EC:4.1.1.19]","[EC:3.5.3.12]","[EC:3.5.1.53]"),]


important2 <- merge(x=important[,1:4],y=filtered_annotated,by.x="product",by.y="product_name",all.x=T)
important <- merge(x=important2[,c(4:9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all.x=T)

#need to do this with grep because names aren't exact
# two<-filtered_annotated[grep('agmatine deiminase',filtered_annotated$product_name,perl=T),]
# three<-filtered_annotated[grep('n-carbamoylputrescine amidase',filtered_annotated$product_name,perl=T),]
# one<-filtered_annotated[grep('arginine decarboxylase',filtered_annotated$product_name,perl=T),]
# four<-rbind(one,two,three)
# important <- merge(x=four[,c(1:5,6,7)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all.x=T)
#important2$genome_name <- lapply(important$genome_name, as.character)
#important[is.na(important)]<-"unknown"

colnames(important)[3:6]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")

#cut out all the unnecessary stuff so they group appropriately
# important$product_name<-gsub('.*arginine decarboxylase.*','arginine decarboxylase',important$product_name,perl = T)
# important$product_name<-gsub('.*agmatine deiminase.*','agmatine deiminase',important$product_name,perl = T)
# important$product_name<-gsub('.*n-carbamoylputrescine amidase.*','n-carbamoylputrescine amidase',important$product_name,perl = T)

important$genome_name<-as.character(important$genome_name)
melted<-melt(important[,c(2:7)])
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("ecnumber","genome","Sample","RNA count")

melted <- with(melted, melted[order(ecnumber, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

#interactive
plot.bar <- ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))
detach(just_poly_products_annot)

#And now for butanoate####
butanoate_annot=read.table("butanoate_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
butanoate_path<-data.frame(unique(butanoate_annot$product))
lowercase_butanoate<-data.frame(tolower(butanoate_path[,1]))
butanoate_path<-lowercase_butanoate
colnames(butanoate_path)<-"V1"
rm(lowercase_butanoate)
just_butanoate_products<-merge(x=sum_by_product_name,y=butanoate_path,by.x="product",by.y="V1",all=F)
butanoate_annot$product <- tolower(butanoate_annot$product)
just_butanoate_products_annot <- merge(x=butanoate_annot,y=just_butanoate_products,by="product",all.x=F,all.y=T)
attach(just_butanoate_products_annot)

#important<-just_butanoate_products_annot[ecnumber %in% c("[EC:2.7.2.7]","[EC:2.3.1.19]"),]
#just buk now
important<-just_butanoate_products_annot[ecnumber %in% c("[EC:2.7.2.7]"),]
important2 <- merge(x=important[,1:4],y=filtered_annotated,by.x="product",by.y="product_name",all.x=T)
#need to do this with grep because names aren't exact
# one<-filtered_annotated[grep('butyrate kinase',filtered_annotated$product_name,perl=T),]
# two<-filtered_annotated[grep('phosphate butyryltransferase',filtered_annotated$product_name,perl=T),]
# three<-rbind(one,two)
important <- merge(x=important2[,c(4:9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all.x=T)
important$genome_name <- lapply(important$genome_name, as.character)
#important[is.na(important)]<-"unknown"


colnames(important)[3:6]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
#important$product_name<-gsub('.*butyrate kinase.*','butyrate kinase',important$product_name,perl = T)
#important$product_name<-gsub('.*phosphate butyryltransferase.*','phosphate butyryltransferase',important$product_name,perl = T)

important$genome_name<-as.character(important$genome_name)
melted<-melt(important[,c(2:7)])
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("ecnumber","genome","Sample","RNA count")

melted <- with(melted, melted[order(ecnumber, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

plot.bar <- ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))

detach(just_butanoate_products_annot)

#for LPS pathway####

lps_annot=read.table("LPS_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
lps_path<-data.frame(unique(lps_annot$product))
lowercase_lps<-data.frame(tolower(lps_path[,1]))
lps_path<-lowercase_lps
colnames(lps_path)<-"V1"
rm(lowercase_lps)
just_lps_products<-merge(x=sum_by_product_name,y=lps_path,by.x="product",by.y="V1",all=F)
lps_annot$product <- tolower(lps_annot$product)
just_lps_products_annot <- merge(x=lps_annot,y=just_lps_products,by="product",all.x=F,all.y=T)
attach(just_lps_products_annot)

important<-just_lps_products_annot[ecnumber %in% c("[EC:3.5.1.108]","[EC:2.3.1.191]"),]
important2 <- merge(x=important[,1:4],y=filtered_annotated,by.x="product",by.y="product_name",all.x=T)
#need to do this with grep because names aren't exact

#lpxC
#one<-filtered_annotated[grep('udp-3-o-\\[3-hydroxymyristoyl\\] n-acetylglucosamine deacetylase',filtered_annotated$ecnumber,perl=T),]

#lpxD
#two<-filtered_annotated[grep('udp-3-o-\\[3-hydroxymyristoyl\\] glucosamine n-acyltransferase',filtered_annotated$ecnumber,perl=T),]
#three<-rbind(one,two)

important <- merge(x=important2[,c(4:9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag",all.x=T)
important$genome_name <- lapply(important$genome_name, as.character)
#important[is.na(important)]<-"unknown"


colnames(important)[3:6]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
#important$ecnumber<-gsub('.*udp-3-o-\\[3-hydroxymyristoyl\\] glucosamine n-acyltransferase.*','UdpOH-Glu-N-Acyl',important$ecnumber,perl = T)
#important$ecnumber<-gsub('.*udp-3-o-\\[3-hydroxymyristoyl\\] n-acetylglucosamine deacetylase.*','UdpOH-N-Acetyl-DeAc',important$ecnumber,perl = T)

important$genome_name<-as.character(important$genome_name)
melted<-melt(important[,c(2:7)])
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("ecnumber","genome","Sample","RNA count")

melted <- with(melted, melted[order(ecnumber, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

plot.bar <- ggplot(data=melted, aes(x=Sample, y=`RNA count`, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~ecnumber) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))

