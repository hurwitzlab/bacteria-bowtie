setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

library(RColorBrewer)
library(reshape2)
library(ggplot2)

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

x=paste(just_lps_products$tracking_id,sep=",",collapse=",")
#copied and pasted x into macvim then %s/,/,\n/g it
#then copied and pasted that into patric feature finder
#at https://www.patricbrc.org/portal/portal/patric/GenomicFeature?cType=taxon&cId=131567&dm=
#and downloaded the results as "FeatureTable.txt"
feature_table = read.delim("FeatureTable.txt",colClasses = "character")
just_lps_products<-merge(x=just_lps_products,y=feature_table[,c("Genome","RefSeq.Locus.Tag")],by.x="tracking_id",by.y="RefSeq.Locus.Tag",all.x=T,all.y=F)
just_lps_products[54,]$Genome<-"Unknown"
write.csv(just_lps_products,"just_lps_products_source_pivot.csv",row.names = F)

melted <- melt(just_lps_products)
jlp_recast <- dcast(melted,melted$product_name + melted$variable ~ melted$Genome,sum)
write.csv(jlp_recast,"jlp_recast.csv",row.names = F)
#[which(just_lps_products$product_name=="3-deoxy-d-manno-octulosonic-acid transferase"

jlp_wide <- reshape(just_lps_products[c("Genome","product_name","S1_FPM","S2_FPM","S3_FPM","S4_FPM","sum")],
                v.names = c("S1_FPM","S2_FPM","S3_FPM","S4_FPM","sum"),
                timevar = "Genome",
                idvar = "product_name",
                direction = "wide")
write.csv(jlp_wide,"jlp_wide.csv")

#this is confusing
#i actually want a clustered, stacked bar chart
#so each cluster will be defined by product_name and S1,S2,S3,S4
#then each stack will be delineated by Genome

#trying again####
jlp<-read.csv("just_lps_products_source_pivot.csv")
jlp<-jlp[,c(2:6,8,9)]

#just lpxC and lpxD
jlp<-read.csv("for first graph.csv")
melted<-melt(jlp)
colnames(melted)<-c("product_name","gene","genome","sample","count")
melted$sample<-gsub("_FPM","",melted$sample)
melted <- with(melted, melted[order(product_name, genome, sample),])
ggplot(data=melted, aes(x=sample, y=count, fill=genome)) + geom_bar(stat="identity") + facet_grid(~gene) + guides(fill=FALSE)
