setwd("/Users/Scott/bacteria-bowtie/scripts/R-interactive")

system(command = "./process-unknown.sh")

setwd("/Users/Scott/unknown-cuffnorm-out")

simple_gene_counts <- read.table("isoforms.fpkm_table",header = T, comment.char = "", strip.white = T, sep = "\t", quote = "", colClasses = c("character","numeric","numeric","numeric","numeric"))

colnames(simple_gene_counts)<-c("tracking_id","S1_FPM","S2_FPM","S3_FPM","S4_FPM")

simple_gene_counts$sum<-rowSums(simple_gene_counts[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])
filtered<-simple_gene_counts[simple_gene_counts$sum!=0,]
rm(simple_gene_counts)

annotation<-read.table("id_to_product.tab",header = F,sep = '\t',quote = "",as.is = T)
lowercase_annotation<-data.frame((annotation[,1]),tolower(annotation[,2]))
annotation<-lowercase_annotation
rm(lowercase_annotation)
colnames(annotation)<-c("tracking_id","product_name")


filtered_annotated<-merge(filtered,annotation,by="tracking_id")

filtered_annotated<-filtered_annotated[order(filtered_annotated$sum,decreasing = T),]

sum_by_product_name<-rowsum(filtered_annotated[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = filtered_annotated$product_name)

sum_by_product_name$product<-row.names(sum_by_product_name)

write.csv(sum_by_product_name,"sum_by_product_name.csv")

write.csv(filtered_annotated,"diff_exp_for_all_bact.csv")
