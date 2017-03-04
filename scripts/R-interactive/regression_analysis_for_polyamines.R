setwd("/Users/Scott/cuffnorm-out")

library(reshape2)
library(tidyr)
library(RColorBrewer)

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

#for polyamines####

polyamines<-read.table("polyamine-search-list",sep="\t",comment.char = "#")
polyamines<-data.frame(tolower(polyamines[,1]))
colnames(polyamines)<-"V1"
polyamines<-unique(polyamines)
poly_product<-data.frame(product=polyamines$V1[1:88])
poly_genes<-data.frame(gene=polyamines$V1[89:203])

just_poly_products<-merge(x=sum_by_product_name,y=poly_product,all=F)
just_poly_genes<-merge(x=sum_by_gene_name,y=poly_genes,all=F)

# mouse expression setup ####
order<-c('s+h+','s-h+')

nos2 <- read.delim("nos2-pcr.tab")

#nos2 <- nos2 %>% mutate(Series =  factor(Series, levels = order)) %>% arrange(Series)

boxplot(nos2$Value~nos2$Series) -> nos2box
nos2box$stats[3,] -> medians_nos2

medians_pcr = data.frame(
  row.names = c("nos2"),
  rbind(medians_nos2)
)
colnames(medians_pcr) = order

just_poly_products[just_poly_products==0]<-NA
polyProductsNoZeroes = na.omit(just_poly_products)

#since we can only correlate on s-h+ and s+h+
colnames(polyProductsNoZeroes) =c('product','s+h-','s-h-','s+h+','s-h+')
row.names(polyProductsNoZeroes)<-polyProductsNoZeroes$product

temp<-polyProductsNoZeroes[,c('s+h+','s-h+')]

polyProductsNoZeroes<-temp
rm(temp)

i=1
for (i in i:(length(medians_pcr[,1]))) {
  j=1
  for (j in j:(length(polyProductsNoZeroes[,1]))) {
    temp1 = as.matrix(medians_pcr[i,])
    temp2 = as.matrix(polyProductsNoZeroes[j,])
    print(paste("Does",row.names(medians_pcr[i,]),"correlate with",row.names(polyProductsNoZeroes[j,]),"?"))
    #print(cor.test(temp1,temp2))
    thing=cor.test(temp1,temp2)
    print(paste("p-value=",thing$p.value))
  }
}

#yup, doesn't work because not enough observations
