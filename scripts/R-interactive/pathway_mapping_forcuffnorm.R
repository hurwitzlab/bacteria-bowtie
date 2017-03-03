#DONE
#these commands make id_to_gene, id_to_product and id_to_ec_number
#setwd("/Users/Scott/bacteria-bowtie/scripts/R-interactive")
#system(command = "./cuffnorm-pathways.sh")
#DONE

library(reshape2)
library(tidyr)
library(RColorBrewer)


setwd("/Users/Scott/cuffnorm-out")

#setup####
filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
filtered_annotated <- filtered_annotated[,2:length(colnames(filtered_annotated))]
attach(filtered_annotated)

#ALREADY DONE
#don't need no hypothetical proteins
# filtered_annotated <- filtered_annotated[grep(".*hypothetical protein.*",product_name,perl=T,invert=T),]
# sum_by_product_name <- read.csv("sum_by_product_name.csv")
# colnames(sum_by_product_name)[1]<-"product"
# attach(sum_by_product_name)

#ALREADY DONE
#don't need no hypothetical proteins
# sum_by_product_name <- sum_by_product_name[grep(".*hypothetical protein.*",product_name,perl=T,invert=T),]
# sum_by_product_name<-sum_by_product_name[,1:5]
# sum_by_product_name$product <- tolower(sum_by_product_name$product)
# sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
# colnames(sum_by_gene_name)[1]<-"gene"
# sum_by_gene_name<-sum_by_gene_name[,1:5]
# sum_by_gene_name$gene <- tolower(sum_by_gene_name$gene)

#now we will attempt to add pathways####
#patric_annotation <- read.delim("all.PATRIC.cds.tab")
#that took over 5 minutes
#probably better to use 'cut' to trim down to the refseq_locus_tag and pathway columns and then load the tab file
# LIKE SO (already done) cut -f 7,21 all.PATRIC.cds.tab > refseq_tag_to_pathway.tab
# and just get lines where pathways are actually known
# grep "\S\t\S" refseq_tag_to_pathway.tab > temp.tab
# mv temp.tab refseq_tag_to_pathway.tab

# Above was for the 1944 set (1944 genomes, not from the year 1944)
# For the combined we have to match product name to pathway colums
# (Already done) cut -f 15,21 all.PATRIC.cds.tab > product_to_pathway.tab
# grep -P '\S\t\S' product_to_pathway.tab > temp.tab
# mv temp.tab product_to_pathway.tab
#patric_annotation <- read.delim("product_to_pathway.tab")
#patric_annotation$product <- tolower(patric_annotation$product)
#with_pathways<-merge(x=filtered_annotated,y=patric_annotation,by.x="product_name",by.y="product")
#This took up 144gb of mem on HPC
#and ended up being a 84gb tab-delimited file!

#ALREADY DONE
#system('cd "/Users/Scott/cuffnorm-out" && cut -f7,19 all.refseq.cds.tab > refseq_locus_to_kegg_pathway.tab')

refseq_annotation <- read.delim("refseq_locus_to_kegg_pathway.tab")
with_pathways<-merge(x=filtered_annotated,y=refseq_annotation,by.x="tracking_id",by.y="refseq_locus_tag",all.x=T)

with_pathways<-separate_rows(with_pathways, pathway, sep = ";")
no_dup<-with_pathways[!duplicated(with_pathways),]
with_pathways<-no_dup
rm(no_dup)

length(with_pathways$pathway[with_pathways$pathway!=''])
#14642 of 146162 actually have annotated pathways, eek!
no_blanks<-with_pathways[with_pathways$pathway!='',]
with_pathways<-no_blanks
rm(no_blanks)

not_na<-with_pathways[!is.na(with_pathways$pathway),]
with_pathways<-not_na
rm(not_na)

write.table(with_pathways,"with_pathways.tab",row.names = F)

sum_by_kegg_pathway<-rowsum(with_pathways[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways$pathway)
sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])

#formatting for stupid bubble plot, which i bet could be done in R
sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)

shortened<-sum_by_kegg_pathway[sum_by_kegg_pathway$sum>mean(sum_by_kegg_pathway$sum),]
shortened<-shortened[order(shortened$sum , decreasing = T),]
shortened<-shortened[,c("Name","S3_FPM","S4_FPM","S1_FPM","S2_FPM")]
colnames(shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")

shortened$combined_effect<-log(shortened$`S-H+ (Combined)`/shortened$`S+H- (Control)`)
shortened$Hhep_effect<-log(shortened$`S+H+ (H. hepaticus only)`/shortened$`S+H- (Control)`)
shortened$smad_effect<-log(shortened$`S-H- (SMAD3 Knockout)`/shortened$`S+H- (Control)`)
shortened<-shortened[order(shortened$combined_effect , decreasing = T),]

#tom's request
shortened$control_effect<-0

effects <- data.matrix(shortened[,6:9])

row.names(effects)<-shortened$Name

x=effects

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

#a heatmap, cuz why not!
heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO","Nothing"))

min(x) #-1.2
max(x) #3.1

new_bubble_source <- shortened[,1:5]

setwd("~/bacteria-bowtie/scripts/R-interactive/")

write.table(new_bubble_source,"sum_by_kegg_pathway_ordered.tab", sep = "\t", quote = T,row.names = F)

system("source ~/.bash_profile && ./bubble.sh sum_by_kegg_pathway_ordered.tab orderedbyeffect")

system("cp orderedbyeffect.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")
