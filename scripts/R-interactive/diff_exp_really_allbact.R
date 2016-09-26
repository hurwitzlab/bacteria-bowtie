setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")

system(command = "./process-reallyall.sh")

setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffdiff_reallyall")

gene_counts <- read.table("total.fpkm_tracking",header = T, comment.char = "", strip.white = T, sep = "\t", quote = "", colClasses = "character")
gene_counts <- gene_counts[order(gene_counts$gene_short_name,decreasing = T),]

simple_gene_counts<-gene_counts[,c("tracking_id","S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")]
simple_gene_counts$sum<-rowSums(simple_gene_counts[,c("S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")])
filtered<-simple_gene_counts[simple_gene_counts$sum!=0,]
filtered<-filtered[order(filtered$sum,decreasing = T),]

annotation<-read.table("id_to_product.tab",header = F,sep = '\t',quote = "")
colnames(annotation)<-c("tracking_id","product_name")

gene_annotation<-read.table("id_to_gene.tab",header = F,sep = '\t',quote = "")
colnames(gene_annotation)<-c("tracking_id","gene")

filtered_annotated<-merge(filtered,annotation,by="tracking_id")
filtered_annotated<-merge(filtered_annotated,gene_annotation,by="tracking_id",all = F)
filtered_annotated<-filtered_annotated[order(filtered_annotated$sum,decreasing = T),]

sum_by_product_name<-rowsum(filtered_annotated[,c("S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")],group = filtered_annotated$product_name)

#write.csv(sum_by_product_name,"sum_by_product_name.csv")

#write.csv(filtered_annotated,"diff_exp_for_all_bact.csv")
