setwd("/Users/Scott/bacteria-bowtie/scripts/R-interactive")

system(command = "./process-combined.sh")

setwd("/Users/Scott/Google Drive/Hurwitz Lab/combined-cuffnorm-out")

simple_gene_counts <- read.table("isoforms.FPKM_table",header = T, comment.char = "", strip.white = T, sep = "\t", quote = "", colClasses = c("character","numeric","numeric","numeric","numeric"))

colnames(simple_gene_counts)<-c("tracking_id","S1_FPM","S2_FPM","S3_FPM","S4_FPM")

simple_gene_counts$sum<-rowSums(simple_gene_counts[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])
filtered<-simple_gene_counts[simple_gene_counts$sum!=0,]

annotation<-read.table("id_to_product.tab",header = F,sep = '\t',quote = "",as.is = T)
lowercase_annotation<-data.frame((annotation[,1]),tolower(annotation[,2]))
annotation<-lowercase_annotation
colnames(annotation)<-c("tracking_id","product_name")
#because parantheses suck and %2c too (,)
gsub('\"','',annotation$product_name)->annotation$product_name
gsub('%2c',',',annotation$product_name)->annotation$product_name

gene_annotation<-read.table("id_to_gene.tab",header = F,sep = '\t',quote = "")
colnames(gene_annotation)<-c("tracking_id","gene")

filtered_annotated<-merge(filtered,annotation,by="tracking_id")
filtered_annotated<-merge(filtered_annotated,gene_annotation,by="tracking_id",all = F)
filtered_annotated<-filtered_annotated[order(filtered_annotated$sum,decreasing = T),]

sum_by_product_name<-rowsum(filtered_annotated[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = filtered_annotated$product_name)
sum_by_gene_name<-rowsum(filtered_annotated[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = filtered_annotated$gene)
sum_by_gene_name$gene<-tolower(row.names(sum_by_gene_name))
sum_by_product_name$product<-row.names(sum_by_product_name)

write.csv(sum_by_product_name,"sum_by_product_name.csv")
write.csv(sum_by_gene_name,"sum_by_gene_name.csv")
write.csv(filtered_annotated,"diff_exp_for_all_bact.csv")
