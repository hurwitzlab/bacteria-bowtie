setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

library(RColorBrewer)

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

#getting more annotation####
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
