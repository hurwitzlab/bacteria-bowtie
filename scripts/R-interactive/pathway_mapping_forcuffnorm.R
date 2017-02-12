#DONE
#these commands make id_to_gene, id_to_product and id_to_ec_number
#setwd("/Users/Scott/tophat-bacteria/scripts/R-interactive")
#system(command = "./cuffnorm-pathways.sh")
#DONE

library(reshape2)
library(tidyr)

setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

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

#now we will attempt to add pathways####
#patric_annotation <- read.delim("all.PATRIC.cds.tab")
#that took over 5 minutes
#probably butter to use 'cut' to trim down to the refseq_locus_tag and pathway columns and then load the tab file
# LIKE SO (already done) cut -f 7,21 all.PATRIC.cds.tab > refseq_tag_to_pathway.tab
# and just get lines where pathways are actually known
# grep "\S\t\S" refseq_tag_to_pathway.tab > temp.tab
# mv temp.tab refseq_tag_to_pathway.tab
patric_annotation <- read.delim("refseq_tag_to_pathway.tab")
with_pathways<-merge(x=filtered_annotated,y=patric_annotation[,c("refseq_locus_tag","pathway")],by.x="tracking_id",by.y="refseq_locus_tag")
with_pathways<-with_pathways[grep(".+",with_pathways$pathway),]
with_pathways<-separate_rows(with_pathways, pathway, sep = ";")

no_dups<-with_pathways[!duplicated(with_pathways),]
with_pathways<-no_dups
rm(no_dups)

sum_by_kegg_pathway<-rowsum(with_pathways[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways$pathway)
sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])

#formatting for stupid bubble plot, which i bet could be done in R
sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)

shortened<-sum_by_kegg_pathway[sum_by_kegg_pathway$sum>mean(sum_by_kegg_pathway$sum),]
shortened<-shortened[order(shortened$sum , decreasing = T),]
shortened<-shortened[,c("Name","S3_FPM","S4_FPM","S1_FPM","S2_FPM")]
colnames(shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")



#setwd("~/bacteria-bowtie/scripts/R-interactive/")
#write.table(shortened,"sum_by_kegg_pathway_above_mean_for_known.tab", sep = "\t", quote = T,row.names = F)

#Have to run this in an external shell cuz ... perl... grumble
#system("./bubble.sh sum_by_kegg_pathway_above_mean.tab Figure7")

#system("cp Figure7.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")

#top_sixty_five<-sum_by_kegg_pathway[order(sum_by_kegg_pathway$sum,decreasing = T)[1:65],]
#top_sixty_five<-top_sixty_five[,c("Name","S1_FPM","S2_FPM","S3_FPM","S4_FPM")]

#write.table(top_sixty_five,"~/tophat-bacteria/scripts/R-interactive/sum_by_kegg_pathway_top_sixty_five.tab", sep = "\t", quote = T,row.names = F)

