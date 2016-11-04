# The purpose of this is to correlate gene expression of the bacteria
# with gene expression of the mouse host
# Specifically genes in the LPS pathway (for bacteria)
# And TLR2/4, NFk-B (for mouse host)
library(RColorBrewer)
library(reshape2)
library(tidyr)

setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

#setup####
filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
filtered_annotated <- filtered_annotated[,2:9]
sum_by_product_name <- read.csv("sum_by_product_name.csv")
colnames(sum_by_product_name)<-c("product","s+h+","s-h+","s+h-","s-h-")
sum_by_product_name<-sum_by_product_name[,1:5]
sum_by_product_name$product <- tolower(sum_by_product_name$product)
sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
colnames(sum_by_gene_name)<-c("gene","s+h+","s-h+","s+h-","s-h-")
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
lps_annot=read.table("LPS_list_annotation",header = T,sep = ";",strip.white = T,quote = "")
lps_annot$product <- tolower(lps_annot$product)
just_lps_products_annot <- merge(x=lps_annot,y=just_lps_products,by="product",all.x=F,all.y=T)

# mouse expression setup ####
tlr4 <- read.delim("tlr4-pcr-data.tab")
order<-c('s+h-','s-h-','s+h+','s-h+')
tlr4$Series <- factor(tlr4$Series,order)
nfkb <- read.delim("nfkb-pcr-data.tab")
nfkb$Series <- factor(nfkb$Series,order)
tlr2 <- read.delim("tlr2-pcr-data.tab")
tlr2$Series <- factor(tlr2$Series,order)

boxplot(tlr4$Values~tlr4$Series) -> tlr4box
boxplot(tlr2$Values~tlr2$Series) -> tlr2box
boxplot(nfkb$Values~nfkb$Series) -> nfkbbox

#these are the medians from boxplot objects
tlr4box$stats[3,] -> medians_tlr4
tlr2box$stats[3,] -> medians_tlr2
nfkbbox$stats[3,] -> medians_nfkb
medians_pcr = data.frame(row.names = c("nfkb","tlr2","tlr4"),rbind(medians_nfkb,medians_tlr2,medians_tlr4))
colnames(medians_pcr) = order

#lets see if these genes correlate with themselves
i=1
for (i in i:(length(medians_pcr[,1])-1)) {
  temp1 = as.matrix(medians_pcr[i,])
  j = i+1
  if (i == 1) {
    for (j in j:(length(medians_pcr[,1]))) {
      temp2 = as.matrix(medians_pcr[j,])
      print(paste(row.names(medians_pcr[i,]),"vs.",row.names(medians_pcr[j,])))
      print(cor.test(temp1,temp2))
    }
  }
  if (i == 2) {
    temp2 = as.matrix(medians_pcr[3,])
    print(paste(row.names(medians_pcr[i,]),"vs.",row.names(medians_pcr[j,])))
    print(cor.test(temp1,temp2))
  }
}

#tlr2 correlates with tlr4

#now trying correlation of mouse pcr to bacterial lps gene products
#first need to remove rows with zeros
just_lps_products[just_lps_products==0]<-NA
lpsProductsNoZeroes = na.omit(just_lps_products)



#old stuff####
boxplot(tlr4$Values ~ tlr4$Series,col='white')

#S-h- vs. S+h-
t.test(tlr4[tlr4$Series=='s+h-',]$Values,tlr4[tlr4$Series=='s-h-',]$Values)

#S-h- vs. S+h+
t.test(tlr4[tlr4$Series=='s-h-',]$Values,tlr4[tlr4$Series=='s+h+',]$Values)

#S+h+ vs. S-h+
t.test(tlr4[tlr4$Series=='s+h+',]$Values,tlr4[tlr4$Series=='s-h+',]$Values)

#S+h- vs. S-h+
t.test(tlr4[tlr4$Series=='s+h-',]$Values,tlr4[tlr4$Series=='s-h+',]$Values)

means_lps <- data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(2997.993999,4566.114029,5776.721804,6506))
means_lps$series <- factor(means_lps$series,order)

#scaling to minimum tlr4 expression
min(tlr4$Values) #0.311
min(means_lps$values) #2997.994
scale_adjust = min(tlr4$Values)-min(means_lps$values)

par(mar = c(5,5,2,5))
boxplot(tlr4$Values ~ tlr4$Series,col='white',ylab='TLR4 relative expression',xlab='Genotype (SMAD3 + or -) and Phenotype (H.hep + or -)')
par(new=T)
plot(x=c(0.5,1.5,2.5,3.5),means_lps$values+scale_adjust, type='l', col='red', axes=F, xlab=NA, ylab=NA,xlim=c(0,4))
axis(side = 4)
mtext(side = 4, line = 3, 'B.fragilis scaled counts')
legend("topleft", legend=c('B.fragilis counts','TLR4 expression'), lty=c(1,0), pch=c(NA, 22), col=c("red", "black"))

by_tlr4<-by(tlr4$Values,tlr4[,"Series"],median)
#by_tlr4<-by(tlr4$Values,tlr4[,"Series"],mean)
medians_tlr4<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_tlr4[[1]],by_tlr4[[2]],by_tlr4[[3]],by_tlr4[[4]]))
cor(medians_tlr4$values,means_lps$values,method = "p")
cor.test(medians_tlr4$values,means_lps$values,method = "p")

b_model <- lm(means_lps$values ~ medians_tlr4$values)
plot(medians_tlr4$values,means_lps$values,ylab='B.fragilis normalized counts',xlab='TLR4 relative expression')
abline(3467,238,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.09862'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))

melted_lps <- melt(just_lps_genes)
by_lps <- by(melted_lps$value,melted_lps[,"variable"],mean)
means_lps<-data.frame(series=c('s+h+','s-h+','s+h-','s-h-'),values=c(by_lps[[1]],by_lps[[2]],by_lps[[3]],by_lps[[4]]))
row.names(means_lps)<-means_lps$series
means_lps<-means_lps[c('s+h-','s-h-','s+h+','s-h+'),]
lin_model<-lm(medians_tlr4$values ~ means_lps$value)
plot(means_lps$values,medians_tlr4$values)
cor(means_lps$values,medians_tlr4$values)
cor.test(means_lps$values,medians_tlr4$values)

boxplot(nfkb$Values ~ nfkb$Series,col='white')

#S-h- vs. S+h-
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s-h-',]$Values)

#S-h- vs. S+h+
t.test(nfkb[nfkb$Series=='s-h-',]$Values,nfkb[nfkb$Series=='s+h+',]$Values)

#S+h+ vs. S-h+
t.test(nfkb[nfkb$Series=='s+h+',]$Values,nfkb[nfkb$Series=='s-h+',]$Values)

#S+h- vs. S-h+
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s-h+',]$Values)

#S+h- vs. S+h+
t.test(nfkb[nfkb$Series=='s+h-',]$Values,nfkb[nfkb$Series=='s+h+',]$Values)

#scaling to min of minimum nfkb expression
min(nfkb$Values) #0.4322686
min(means_lps$values) #2997.994
scale_adjust = min(nfkb$Values)-min(means_lps$values)
par(mar = c(5,5,2,5))
boxplot(nfkb$Values ~ nfkb$Series,col='white',ylab='NFkB relative expression',xlab='Genotype (SMAD3 + or -) and Phenotype (H.hep + or -)')
par(new=T)
plot(x=c(0.5,1.5,2.5,3.5),means_lps$values+scale_adjust, type='l', col='red', axes=F, xlab=NA, ylab=NA,xlim=c(0,4))
axis(side = 4)
mtext(side = 4, line = 3, 'B.fragilis scaled counts')
legend("topleft", legend=c('B.fragilis counts','NFkB expression'), lty=c(1,0), pch=c(NA, 22), col=c("red", "black"))


clusterof2 <- rbind(nfkb,tlr4)
clusterof2$Series <- factor(clusterof2$Series,order)
boxplot(Values ~ Series, data = clusterof2)

by_clusterof2<-by(clusterof2$Values,clusterof2[,"Series"],median)
#by_clusterof2<-by(clusterof2$Values,clusterof2[,"Series"],mean)
medians_clusterof2<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_clusterof2[[1]],by_clusterof2[[2]],by_clusterof2[[3]],by_clusterof2[[4]]))
cor(medians_clusterof2$values,means_lps$values,method = "p")
cor.test(medians_clusterof2$values,means_lps$values,method = "p")

b_model <- lm(means_lps$values ~ medians_clusterof2$values)
plot(medians_clusterof2$values,means_lps$values,ylab='B.fragilis normalized counts',xlab='TLR4+nfkb median expression')
abline(1620,1323,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.01459'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))

boxplot(tlr2$Values ~ tlr2$Series,col='white')

clusterof3 <- rbind(nfkb,tlr4,tlr2)
clusterof3$Series <- factor(clusterof3$Series,order)
boxplot(Values ~ Series, data = clusterof3)

by_clusterof3<-by(clusterof3$Values,clusterof3[,"Series"],median)
#by_clusterof3<-by(clusterof3$Values,clusterof3[,"Series"],mean)
medians_clusterof3<-data.frame(series=c('s+h-','s-h-','s+h+','s-h+'),values=c(by_clusterof3[[1]],by_clusterof3[[2]],by_clusterof3[[3]],by_clusterof3[[4]]))
cor(medians_clusterof3$values,means_lps$values,method = "p")
cor.test(medians_clusterof3$values,means_lps$values,method = "p")

b_model <- lm(means_lps$values ~ medians_clusterof3$values)
plot(medians_clusterof3$values,means_lps$values,ylab='B.fragilis normalized counts',xlab='TLR4/TLR2/nfkb relative expression')
abline(2205,1128,col='blue')
legend("bottomright", legend=c('Line of regression','Correlated values', 'pvalue = 0.1247'), lty=c(1,0,0), pch=c(NA, 1,-1), col=c("blue", "black", "black"))
