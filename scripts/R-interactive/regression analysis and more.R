# The purpose of this is to correlate gene expression of the bacteria
# with gene expression of the mouse host
# Specifically genes in the LPS pathway (for bacteria)
# And TLR2/4, NFk-B (for mouse host)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(dplyr)

setwd("/Users/Scott/combined-cuffnorm-out")

#setup####
# filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
# filtered_annotated <- filtered_annotated[,2:9]
# sum_by_product_name <- read.csv("sum_by_product_name.csv")
# colnames(sum_by_product_name)<-c("product","s+h+","s-h+","s+h-","s-h-")
# sum_by_product_name<-sum_by_product_name[,1:5]
# sum_by_product_name$product <- tolower(sum_by_product_name$product)
# sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
# colnames(sum_by_gene_name)<-c("gene","s+h+","s-h+","s+h-","s-h-")
# sum_by_gene_name<-sum_by_gene_name[,1:5]
# sum_by_gene_name$gene <- tolower(sum_by_gene_name$gene)

#lucky for us, we now have the full LPS pathway already parsed out via the "pathway_mapping_for_combined.R"
lps<-read.table("all_lps_products.tab",header=T)
lps_nodup<-lps[!duplicated(lps),]
rm(lps)
# lpxC=lps_nodup[grep('.*3.5.1.108.*',lps_nodup$product_name,perl=T),]
# lpxC$gene='lpxC'
sum_by_product_name<-rowsum(lps_nodup[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = lps_nodup$product_name)


#for LPS pathway####
# lps_path<-as.character(lps_nodup$product_name[!duplicated(lps_nodup)])
# lowercase_lps<-data.frame(tolower(lps_path[,1]))
# lps_path<-lowercase_lps
# colnames(lps_path)<-"V1"
# rm(lowercase_lps)
#lps_genes<-data.frame(gene=lps_path[34:96,])
# lps_products<-data.frame(product=lps_path[,])
#
# just_lps_products<-merge(x=sum_by_product_name,y=lps_products,by="product",all=F)
# just_lps_genes<-merge(x=sum_by_gene_name,y=lps_genes,by="gene",all=F)
# lps_annot=read.table("LPS_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
# lps_annot$product <- tolower(lps_annot$product)
# just_lps_products_annot <- merge(x=lps_annot,y=just_lps_products,by="product",all.x=F,all.y=T)


# mouse expression setup ####
order<-c('s+h-','s-h-','s+h+','s-h+')

tlr4 <- read.delim("tlr4-pcr-data.tab")
tlr4 <- tlr4 %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)
nfkb <- read.delim("nfkb-pcr-data.tab")
nfkb <- nfkb %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)
tlr2 <- read.delim("tlr2-pcr-data.tab")
tlr2 <- tlr2 %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)
#another inflammatory pathway gene in mice
irak4 <- read.delim("irak4-pcr-data.tab")
irak4 <- irak4 %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)
#and another
cd14 <- read.delim("CD14-pcr-data.tab")
cd14 <- cd14 %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)

boxplot(tlr4$Values~tlr4$Series) -> tlr4box
boxplot(tlr2$Values~tlr2$Series) -> tlr2box
boxplot(nfkb$Values~nfkb$Series) -> nfkbbox
boxplot(irak4$Values~irak4$Series) -> irak4box
boxplot(cd14$Values~cd14$Series) -> cd14box

#these are the medians from boxplot objects
tlr4box$stats[3,] -> medians_tlr4
tlr2box$stats[3,] -> medians_tlr2
nfkbbox$stats[3,] -> medians_nfkb
irak4box$stats[3,] -> medians_irak4
cd14box$stats[3,] -> medians_cd14
medians_pcr = data.frame(
  row.names = c("nfkb","tlr2","tlr4","irak4","cd14"),
  rbind(medians_nfkb,medians_tlr2,medians_tlr4,medians_irak4,
        medians_cd14)
  )
colnames(medians_pcr) = order

#lets see if these genes correlate with themselves
# i=1
# for (i in i:(length(medians_pcr[,1])-1)) {
#   temp1 = as.matrix(medians_pcr[i,])
#   j = i+1
#   if (i == 1) {
#     for (j in j:(length(medians_pcr[,1]))) {
#       temp2 = as.matrix(medians_pcr[j,])
#       print(paste(row.names(medians_pcr[i,]),"vs.",row.names(medians_pcr[j,])))
#       print(cor.test(temp1,temp2))
#     }
#   }
#   if (i == 2) {
#     temp2 = as.matrix(medians_pcr[3,])
#     print(paste(row.names(medians_pcr[i,]),"vs.",row.names(medians_pcr[j,])))
#     print(cor.test(temp1,temp2))
#   }
# }

#tlr2 correlates with tlr4

#now trying correlation of mouse pcr to bacterial lps gene products
#first need to remove rows with zeros
sum_by_product_name[sum_by_product_name==0]<-NA
lpsProductsNoZeroes = na.omit(sum_by_product_name)

colnames(lpsProductsNoZeroes) =c('s+h-','s-h-','s+h+','s-h+')

lpxC=lpsProductsNoZeroes[grep('.*3.5.1.108.*',row.names(lpsProductsNoZeroes),perl=T),]
lpxD=lpsProductsNoZeroes[grep('.*2.3.1.191.*',row.names(lpsProductsNoZeroes),perl=T),]
lpxK=lpsProductsNoZeroes[grep('.*2.7.1.130.*',row.names(lpsProductsNoZeroes),perl=T),]
kdtA=lpsProductsNoZeroes[grep('.*2.4.99.*',row.names(lpsProductsNoZeroes),perl=T),]

#row.names(lpsProductsNoZeroes) = c("kdtA, waaA","waaL/rfaL","lpxD","lpxC")

i=1
for (i in i:(length(medians_pcr[,1]))) {
  j=1
  for (j in j:(length(lpsProductsNoZeroes[,1]))) {
    temp1 = as.matrix(medians_pcr[i,])
    temp2 = as.matrix(lpsProductsNoZeroes[j,])
    print(paste("Does",row.names(medians_pcr[i,]),"correlate with",row.names(lpsProductsNoZeroes[j,]),"?"))
    #print(cor.test(temp1,temp2))
    thing=cor.test(temp1,temp2)
    print(paste("p-value=",thing$p.value))
  }
}

#Interesting, in this new set
#TLR2 correlates with kdtA/waaA and lpxD
#While TLR4 correlates with lpxD and lpxC
#Have no idea if this really means anything
#They all are close to significance

#correlation plots####
#in order of p-value results

#kdtA/waaA and tlr2
plot(unlist(medians_pcr["tlr2",]),unlist(lpsProductsNoZeroes["kdtA, waaA",]),xlab="Relative expression of TLR2",ylab="Relative expression of kdtA/waaA",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.01726")
model<-lm(unlist(lpsProductsNoZeroes["kdtA, waaA",]) ~ unlist(medians_pcr["tlr2",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)

#lpxD and tlr2
plot(unlist(medians_pcr["tlr2",]),unlist(lpsProductsNoZeroes["lpxD",]),xlab="Relative expression of TLR2",ylab="Relative expression of lpxD",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.026")

model<-lm(unlist(lpsProductsNoZeroes["lpxD",]) ~ unlist(medians_pcr["tlr2",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)

#lpxC and tlr2
plot(unlist(medians_pcr["tlr2",]),unlist(lpsProductsNoZeroes["lpxC",]),xlab="Relative expression of TLR2",ylab="Relative expression of lpxC",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.06263")

model<-lm(unlist(lpsProductsNoZeroes["lpxC",]) ~ unlist(medians_pcr["tlr2",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)

#kdtA/waaA and tlr4
plot(unlist(medians_pcr["tlr4",]),unlist(lpsProductsNoZeroes["kdtA, waaA",]),xlab="Relative expression of TLR4",ylab="Relative expression of kdtA/waaA",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.06248")

model<-lm(unlist(lpsProductsNoZeroes["kdtA, waaA",]) ~ unlist(medians_pcr["tlr4",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)

#lpxD and tlr4
plot(unlist(medians_pcr["tlr4",]),unlist(lpsProductsNoZeroes["lpxD",]),xlab="Relative expression of TLR4",ylab="Relative expression of lpxD",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.004623")
model<-lm(unlist(lpsProductsNoZeroes["lpxD",]) ~ unlist(medians_pcr["tlr4",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)

#lpxC and tlr4
plot(unlist(medians_pcr["tlr4",]),unlist(lpsProductsNoZeroes["lpxC",]),xlab="Relative expression of TLR4",ylab="Relative expression of lpxC",cex.axis=1.5,cex.lab=1.5,pch=21,col="black",bg="black",main="p-value=0.03969")
model<-lm(unlist(lpsProductsNoZeroes["lpxC",]) ~ unlist(medians_pcr["tlr4",]))
abline(model$coefficients[[1]],model$coefficients[[2]],col="blue",lwd=2)


