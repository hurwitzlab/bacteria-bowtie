setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

library(reshape2)
library(tidyr)
library(RColorBrewer)

#setup####
with_pathways<-read.table("with_pathways.tab",header=T)

#oxidative phosphorylation####

attach(with_pathways)

oxphos<-with_pathways[pathway=="00190|Oxidative phosphorylation",]
no_blanks<-oxphos[oxphos$product_name!="",]
oxphos<-no_blanks
rm(no_blanks)

oxphos_product_sums<-rowsum(oxphos[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group=oxphos$product_name)
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

#for LPS pathway####
lps_path<-read.table("LPS_search_list",sep = "\t",comment.char = "#")
lowercase_lps<-data.frame(tolower(lps_path[,1]))
lps_path<-lowercase_lps
colnames(lps_path)<-"V1"
rm(lowercase_lps)
lps_genes<-data.frame(gene=lps_path[34:96,])
lps_products<-data.frame(product=lps_path[1:34,])

just_lps_products<-merge(x=sum_by_product_name,y=lps_products,by="product",all=F)
just_lps_genes<-merge(x=sum_by_gene_name,y=lps_genes,by="gene",all=F)

#for polyamines####
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
# write.csv(just_poly_products_annot,"poly_products_annotated.csv")
#more diffexp stuff / heatmap####

just_poly_from_excel <- read.csv("poly_products_for_kegg_figure.csv",header=T)

row.names(just_poly_from_excel)<-just_poly_from_excel$id_on_kegg

just_poly_from_excel <- data.matrix(just_poly_from_excel[,3:5])

x=just_poly_from_excel

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none",margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))

# and now heatmap for LPS####

lps_annot=read.table("LPS_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
lps_annot$product <- tolower(lps_annot$product)
just_lps_products_annot <- merge(x=lps_annot,y=just_lps_products,by="product",all.x=F,all.y=T)
#write.csv(just_lps_products_annot,"just_lps_products_annot.csv")

just_lps_from_excel <- read.delim("for_lps_heatmap.txt",header = T)

row.names(just_lps_from_excel)<-just_lps_from_excel$ecnumber

just_lps_from_excel <- just_lps_from_excel[order(just_lps_from_excel$cancer),]

just_lps_from_excel <- data.matrix(just_lps_from_excel[,4:6])

x=just_lps_from_excel

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))

just_lps_from_excel <- read.delim("for_lps_heatmap.txt",header = T)

row.names(just_lps_from_excel)<-just_lps_from_excel$gene

just_lps_from_excel <- just_lps_from_excel[order(just_lps_from_excel$cancer),]

just_lps_from_excel <- data.matrix(just_lps_from_excel[,4:6])

x=just_lps_from_excel

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))

#for Butanoate pathway####
butanoate_annot=read.table("butanoate_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
butanoate_path<-data.frame(unique(butanoate_annot$product))
lowercase_butanoate<-data.frame(tolower(butanoate_path[,1]))
butanoate_path<-lowercase_butanoate
colnames(butanoate_path)<-"V1"
rm(lowercase_butanoate)
# butanoate_genes<-data.frame(gene=butanoate_path[34:96,])
# butanoate_products<-data.frame(product=butanoate_path[1:34,])

just_butanoate_products<-merge(x=sum_by_product_name,y=butanoate_path,by.x="product",by.y="V1",all=F)
# just_butanoate_genes<-merge(x=sum_by_gene_name,y=butanoate_genes,by="gene",all=F)


# and now heatmap for butanoate####


butanoate_annot$product <- tolower(butanoate_annot$product)
just_butanoate_products_annot <- merge(x=butanoate_annot,y=just_butanoate_products,by="product",all.x=F,all.y=T)
#write.csv(just_butanoate_products_annot,"just_butanoate_products_annot.csv")

just_butanoate_from_excel <- read.delim("for_butanoate_heatmap.txt",header = T)

row.names(just_butanoate_from_excel)<-just_butanoate_from_excel$ecnumber

just_butanoate_from_excel <- just_butanoate_from_excel[order(just_butanoate_from_excel$Combined),]

just_butanoate_from_excel <- data.matrix(just_butanoate_from_excel[,4:6])

x=just_butanoate_from_excel

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))

just_butanoate_from_excel <- read.delim("for_butanoate_heatmap.txt",header = T)

row.names(just_butanoate_from_excel)<-just_butanoate_from_excel$gene

just_butanoate_from_excel <- just_butanoate_from_excel[order(just_butanoate_from_excel$Combined),]

just_butanoate_from_excel <- data.matrix(just_butanoate_from_excel[,4:6])

x=just_butanoate_from_excel

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))
