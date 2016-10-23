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

#for m schaedleri####
mschaedleri <- read.csv("MschaedleriFeatures.csv")
mschaedleri_exp <- merge(mschaedleri,filtered_annotated,by.x="RefSeq.Locus.Tag",by.y="tracking_id",all=F)
mschaedleri_exp_vs_all_products <- merge(mschaedleri_exp,sum_by_product_name,by.x="product_name",by.y="product",all=F)

#write.csv(mschaedleri_exp_vs_all_products,"mschaedleri_contribution_to_set.csv")


