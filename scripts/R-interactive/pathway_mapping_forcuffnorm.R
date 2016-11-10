#DONE
#these commands make id_to_gene, id_to_product and id_to_ec_number
#setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")
#system(command = "./cuffnorm-pathways.sh")
#DONE

library(reshape2)
library(tidyr)

setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

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

#now we will attempt to add pathways####
#patric_annotation <- read.delim("all.PATRIC.cds.tab")
#that took over 5 minutes
#probably butter to use 'cut' to trim down to the refseq_locus_tag and pathway columns and then load the tab file
# LIKE SO (already done) cut -f 7,21 all.PATRIC.cds.tab > refseq_tag_to_pathway.tab
# and just get lines where pathways are actually known
# grep "\S\t\S" refseq_tag_to_pathway.tab > temp.tab
# mv temp.tab refseq_tag_to_pathway.tab
patric_annotation <- read.delim("refseq_tag_to_pathway.tab")
with_pathways<-merge(x=filtered_annotated,y=patric_annotation[,c("refseq_locus_tag","pathway")],by.x="tracking_id",by.y="refseq_locus_tag")
with_pathways<-with_pathways[grep(".+",with_pathways$pathway),]
with_pathways<-separate_rows(with_pathways, pathway, sep = ";")
sum_by_kegg_pathway<-rowsum(with_pathways[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways$pathway)
sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])

#formatting for stupid bubble plot, which i bet could be done in R
sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)

shortened<-sum_by_kegg_pathway[sum_by_kegg_pathway$sum>mean(sum_by_kegg_pathway$sum),]
shortened<-shortened[order(shortened$sum , decreasing = T),]
shortened<-shortened[,c("Name","S1_FPM","S2_FPM","S3_FPM","S4_FPM")]
top_sixty_five<-sum_by_kegg_pathway[order(sum_by_kegg_pathway$sum,decreasing = T)[1:65],]
top_sixty_five<-top_sixty_five[,c("Name","S1_FPM","S2_FPM","S3_FPM","S4_FPM")]

write.table(top_sixty_five,"~/tophat-bacteria/scripts/R-interactive/sum_by_kegg_pathway_top_sixty_five.tab", sep = "\t", quote = T,row.names = F)

#setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")

#write.csv(top_sixty_five,"sum_by_kegg_pathway_top_sixty_five.csv",quote = T,row.names = F)

# Maybe use later
# #for LPS pathway####
# lps_path<-read.table("LPS_search_list",sep = "\t",comment.char = "#")
# lowercase_lps<-data.frame(tolower(lps_path[,1]))
# lps_path<-lowercase_lps
# colnames(lps_path)<-"V1"
# rm(lowercase_lps)
# lps_genes<-data.frame(gene=lps_path[34:96,])
# lps_products<-data.frame(product=lps_path[1:34,])
#
# just_lps_products<-merge(x=sum_by_product_name,y=lps_products,by="product",all=F)
# just_lps_genes<-merge(x=sum_by_gene_name,y=lps_genes,by="gene",all=F)
#
# #for polyamines####
# polyamines<-read.table("polyamine-search-list",sep="\t",comment.char = "#")
# polyamines<-data.frame(tolower(polyamines[,1]))
# colnames(polyamines)<-"V1"
# polyamines<-unique(polyamines)
# poly_product<-data.frame(product=polyamines$V1[1:88])
# poly_genes<-data.frame(gene=polyamines$V1[89:203])
#
# just_poly_products<-merge(x=sum_by_product_name,y=poly_product,all=F)
# just_poly_genes<-merge(x=sum_by_gene_name,y=poly_genes,all=F)
#
# #getting more annotation for poly####
# poly_annot<-read.table("polyamine_list_annotation",header = T,sep = ";",strip.white = T)
# poly_annot$product <- tolower(poly_annot$product)
# just_poly_products_annot <- merge(x=poly_annot,y=just_poly_products)
# # write.csv(just_poly_products_annot,"poly_products_annotated.csv")
#
# #more diffexp stuff / heatmap####
#
# just_poly_from_excel <- read.csv("poly_products_for_kegg_figure.csv",header=T)
#
# row.names(just_poly_from_excel)<-just_poly_from_excel$id_on_kegg
#
# just_poly_from_excel <- data.matrix(just_poly_from_excel[,3:5])
#
# x=just_poly_from_excel
#
# oldPar <- par(no.readonly = T)
#
# myColors=colorRampPalette(c("Blue","Yellow"))
#
# heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none",margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))
#
# # and now heatmap for LPS####
#
# lps_annot=read.table("LPS_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
# lps_annot$product <- tolower(lps_annot$product)
# just_lps_products_annot <- merge(x=lps_annot,y=just_lps_products,by="product",all.x=F,all.y=T)
# #write.csv(just_lps_products_annot,"just_lps_products_annot.csv")
#
# just_lps_from_excel <- read.delim("for_lps_heatmap.txt",header = T)
#
# row.names(just_lps_from_excel)<-just_lps_from_excel$ecnumber
#
# just_lps_from_excel <- just_lps_from_excel[order(just_lps_from_excel$cancer),]
#
# just_lps_from_excel <- data.matrix(just_lps_from_excel[,4:6])
#
# x=just_lps_from_excel
#
# oldPar <- par(no.readonly = T)
#
# myColors=colorRampPalette(c("Blue","Yellow"))
#
# heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))
#
# just_lps_from_excel <- read.delim("for_lps_heatmap.txt",header = T)
#
# row.names(just_lps_from_excel)<-just_lps_from_excel$gene
#
# just_lps_from_excel <- just_lps_from_excel[order(just_lps_from_excel$cancer),]
#
# just_lps_from_excel <- data.matrix(just_lps_from_excel[,4:6])
#
# x=just_lps_from_excel
#
# oldPar <- par(no.readonly = T)
#
# myColors=colorRampPalette(c("Blue","Yellow"))
#
# heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Cancer", "Inflammation","SMAD3-KO"))
#
# #for Butanoate pathway####
# butanoate_annot=read.table("butanoate_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
# butanoate_path<-data.frame(unique(butanoate_annot$product))
# lowercase_butanoate<-data.frame(tolower(butanoate_path[,1]))
# butanoate_path<-lowercase_butanoate
# colnames(butanoate_path)<-"V1"
# rm(lowercase_butanoate)
# # butanoate_genes<-data.frame(gene=butanoate_path[34:96,])
# # butanoate_products<-data.frame(product=butanoate_path[1:34,])
#
# just_butanoate_products<-merge(x=sum_by_product_name,y=butanoate_path,by.x="product",by.y="V1",all=F)
# # just_butanoate_genes<-merge(x=sum_by_gene_name,y=butanoate_genes,by="gene",all=F)
#
# # and now heatmap for butanoate####
#
#
# butanoate_annot$product <- tolower(butanoate_annot$product)
# just_butanoate_products_annot <- merge(x=butanoate_annot,y=just_butanoate_products,by="product",all.x=F,all.y=T)
# write.csv(just_butanoate_products_annot,"just_butanoate_products_annot.csv")
#
# just_butanoate_from_excel <- read.delim("for_butanoate_heatmap.txt",header = T)
#
# row.names(just_butanoate_from_excel)<-just_butanoate_from_excel$ecnumber
#
# just_butanoate_from_excel <- just_butanoate_from_excel[order(just_butanoate_from_excel$Combined),]
#
# just_butanoate_from_excel <- data.matrix(just_butanoate_from_excel[,4:6])
#
# x=just_butanoate_from_excel
#
# oldPar <- par(no.readonly = T)
#
# myColors=colorRampPalette(c("Blue","Yellow"))
#
# heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))
#
# just_butanoate_from_excel <- read.delim("for_butanoate_heatmap.txt",header = T)
#
# row.names(just_butanoate_from_excel)<-just_butanoate_from_excel$gene
#
# just_butanoate_from_excel <- just_butanoate_from_excel[order(just_butanoate_from_excel$Combined),]
#
# just_butanoate_from_excel <- data.matrix(just_butanoate_from_excel[,4:6])
#
# x=just_butanoate_from_excel
#
# oldPar <- par(no.readonly = T)
#
# myColors=colorRampPalette(c("Blue","Yellow"))
#
# heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))
