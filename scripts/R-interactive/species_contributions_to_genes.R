setwd("/Users/Scott/Google Drive/Hurwitz Lab/combined-cuffnorm-out")

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

#for m schaedleri####
# mschaedleri <- read.csv("MschaedleriFeatures.csv")
# mschaedleri_exp <- merge(mschaedleri,filtered_annotated,by.x="RefSeq.Locus.Tag",by.y="tracking_id",all=F)
# mschaedleri_exp_vs_all_products <- merge(mschaedleri_exp,sum_by_product_name,by.x="product_name",by.y="product",all=F)

#write.csv(mschaedleri_exp_vs_all_products,"mschaedleri_contribution_to_set.csv")

#for LPS pathway####
lps_path<-read.table("LPS_search_list",sep = "\t",comment.char = "#")
lps_genes<-data.frame(gene=lps_path[34:96,])
lowercase_lps<-data.frame(tolower(lps_path[,1]))
lps_path<-lowercase_lps
colnames(lps_path)<-"V1"
rm(lowercase_lps)
lps_products<-data.frame(product=lps_path[1:34,])

just_lps_products<-merge(x=filtered_annotated,y=lps_products,by.x="product_name",by.y="product",all=F)
just_lps_genes<-merge(x=filtered_annotated,y=lps_genes,by="gene",all=F)

# x=paste(just_lps_products$tracking_id,sep=",",collapse=",")
# #copied and pasted x into macvim then %s/,/,\n/g it
# #then copied and pasted that into patric feature finder
# #at https://www.patricbrc.org/portal/portal/patric/GenomicFeature?cType=taxon&cId=131567&dm=
# #and downloaded the results as "FeatureTable.txt"
# feature_table = read.delim("FeatureTable.txt",colClasses = "character")
# just_lps_products<-merge(x=just_lps_products,y=feature_table[,c("Genome","RefSeq.Locus.Tag")],by.x="tracking_id",by.y="RefSeq.Locus.Tag",all.x=T,all.y=F)
# just_lps_products[54,]$Genome<-"Unknown"
# write.csv(just_lps_products,"just_lps_products_source_pivot.csv",row.names = F)
#
# melted <- melt(just_lps_products)
# jlp_recast <- dcast(melted,melted$product_name + melted$variable ~ melted$Genome,sum)
# write.csv(jlp_recast,"jlp_recast.csv",row.names = F)
# #[which(just_lps_products$product_name=="3-deoxy-d-manno-octulosonic-acid transferase"
#
# jlp_wide <- reshape(just_lps_products[c("Genome","product_name","S1_FPM","S2_FPM","S3_FPM","S4_FPM","sum")],
#                 v.names = c("S1_FPM","S2_FPM","S3_FPM","S4_FPM","sum"),
#                 timevar = "Genome",
#                 idvar = "product_name",
#                 direction = "wide")
# write.csv(jlp_wide,"jlp_wide.csv")

#this is confusing
#i actually want a clustered, stacked bar chart
#so each cluster will be defined by product_name and S1,S2,S3,S4
#then each stack will be delineated by Genome

#trying again####
#jlp<-read.csv("just_lps_products_source_pivot.csv")
#jlp<-jlp[,c(2:6,8,9)]

#just lpxC and lpxD
jlp<-read.csv("for first graph.csv")
colnames(jlp)[2:5]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
melted<-melt(jlp)
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("product_name","gene","genome","Sample","count")
melted <- with(melted, melted[order(product_name, genome, Sample),])

#Now interactive!!!!
plot.bar <- ggplot(data=melted, aes(x=Sample, y=count, fill=genome, tooltip = genome, data_id = Sample)) + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

ggiraph(code = print(plot.bar))

plot.bar + geom_bar_interactive(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)

#just kdtA and waaL
jlp<-read.csv("for second graph.csv")
colnames(jlp)[2:5]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
melted<-melt(jlp)
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("product_name","gene","genome","Sample","count")
melted <- with(melted, melted[order(product_name, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=count, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)

#Let's try polyamines now####
#Got all.PATRIC.cds.tab from the server and cut down to just genome names and refseq locus tags
#system("cut -f2,7 all.PATRIC.cds.tab > genome-name_to_refseq-locus-tag")
genome_to_feature <- read.delim("genome-name_to_refseq-locus-tag")

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
important2 <- merge(x=important[,1:4],y=filtered_annotated,by.x="product",by.y="product_name",all.x=T)
important <- merge(x=important2[,c(1,3,5:9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag")
colnames(important)[4:7]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
melted<-melt(important[,2:8])
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("product_name","gene","genome","Sample","count")

melted <- with(melted, melted[order(product_name, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=count, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)
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

important<-just_butanoate_products_annot[ecnumber %in% c("[EC:2.7.2.7]","[EC:2.3.1.19]"),]
important2 <- merge(x=important[,1:4],y=filtered_annotated,by.x="product",by.y="product_name",all.x=T)
important <- merge(x=important2[,c(1,3,5:9)],y=genome_to_feature,by.x="tracking_id",by.y="refseq_locus_tag")
colnames(important)[4:7]<-c("S+H+ (H. hepaticus only)","S-H+ (Combined)","S+H- (Control)","S-H- (SMAD3 Knockout)")
melted<-melt(important[,2:8])
order<-c("S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
melted <- melted %>% mutate(variable =  factor(variable, levels = order)) %>% arrange(variable)
colnames(melted)<-c("product_name","gene","genome","Sample","count")

melted <- with(melted, melted[order(product_name, genome, Sample),])
plot.bar = ggplot(data=melted, aes(x=Sample, y=count, fill=genome))
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ guides(fill=FALSE)
plot.bar + geom_bar(stat="identity", col="black", size = .5) + facet_grid(~gene) + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)
