setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")

system(command = "./process-mouse.sh")

setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffdiff_mouse_v1/")
#
# library(cummeRbund)
#
# sqliteQuickSQL<-dbGetQuery
#
# dbBeginTransaction<-dbBegin
#
#cuff=readCufflinks(dbFile="cuffData.db",genome="multiple-bacteria",reload=T)

gene_diff <- read.table("isoform_exp.diff",header = T)
gene_diff <- gene_diff[order(gene_diff$gene,decreasing = T),]
gene_counts <- read.table("isoforms.fpkm_tracking",header = T)
gene_counts <- gene_counts[order(gene_counts$gene_short_name,decreasing = T),]

simple_gene_counts<-gene_counts[,c("tracking_id","X.S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")]
colnames(simple_gene_counts)<-c("tracking_id","S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")
simple_gene_counts$sum<-rowSums(simple_gene_counts[,c("S1_FPKM","S2_FPKM","S3_FPKM","S4_FPKM")])
filtered<-simple_gene_counts[simple_gene_counts$sum!=0,]
filtered<-filtered[order(filtered$sum,decreasing = T),]

gene_annotation<-read.table("id_to_gene.tab",header = F,sep = '\t',quote = "")
colnames(gene_annotation)<-c("tracking_id","gene")

filtered_annotated<-merge(filtered,gene_annotation,by="tracking_id",all = F)
filtered_annotated<-filtered_annotated[order(filtered_annotated$sum,decreasing = T),]

write.csv(filtered_annotated,"diff_exp_for_mouse.csv")
