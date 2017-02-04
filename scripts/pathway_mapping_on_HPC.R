#This should be run using Rscript and takes a command-line
#argument of a directory

library(reshape2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

setwd(paste(args[1]))

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
#because parantheses suck and %2c too (,)
gsub('\"','',annotation$product_name)->annotation$product_name
gsub('%2c',',',annotation$product_name)->annotation$product_name
#trying to remove some of the inconsistent names from products
#and hypotheticals, probables, predicted
annotation<-annotation[grep(".*hypothetical protein.*",annotation$product_name,perl=T,invert=T),]
annotation<-annotation[grep(".*probable.*",annotation$product_name,perl=T,invert=T),]
annotation<-annotation[grep(".*predicted.*",annotation$product_name,perl=T,invert=T),]
annotation<-annotation[grep(".*uncharacterized.*",annotation$product_name,perl=T,invert=T),]

annotation<-annotation[grep(".*putative*",annotation$product_name,perl=T,invert=T),]

known_species<-annotation[grep('^fig\\|6666666.*',annotation$tracking_id,perl=T,invert=T),]
unknown_species<-annotation[grep('^fig\\|6666666.*',annotation$tracking_id,perl=T),]
print("making annotation_best_i_can_do")

#this filters out product names that are not consistent with known species
matrix_of_goodness <- merge(known_species,unknown_species,by="product_name",all.x=T)
good_unknown <- unknown_species[unknown_species$tracking_id %in% matrix_of_goodness$tracking_id.y,]
annotation_best_i_can_do<-rbind(good_unknown,known_species)

write.table(x = good_unknown,"good_gene_names_from_unknown.tab")
#
# rm(annotation,known_species,unknown_species,matrix_of_goodness)
#
# #f$%^ genes
# #gene_annotation<-read.table("id_to_gene.tab",header = F,sep = '\t',quote = "")
# #colnames(gene_annotation)<-c("tracking_id","gene")
#
# print("making filtered_annotated")
# #print(colnames(filtered))
# #print(colnames(annotation_best_i_can_do))
# filtered_annotated<-merge(filtered,annotation_best_i_can_do,by="tracking_id")
# filtered_annotated<-merge(filtered_annotated,gene_annotation,by="tracking_id",all.x=T,sort=F)
# filtered_annotated<-filtered_annotated[order(filtered_annotated$sum,decreasing = T),]
#
# sum_by_product_name<-rowsum(filtered_annotated[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = filtered_annotated$product_name)
# sum_by_gene_name<-rowsum(filtered_annotated[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = filtered_annotated$gene)
# sum_by_gene_name$gene<-tolower(row.names(sum_by_gene_name))
# sum_by_product_name$product<-row.names(sum_by_product_name)
#
# #don't need row names as a column
# write.csv(sum_by_product_name,"sum_by_product_name.csv",row.names = F)
# write.csv(sum_by_gene_name,"sum_by_gene_name.csv",row.names = F)
# write.csv(filtered_annotated,"diff_exp_for_all_bact.csv",row.names = F)
# rm(sum_by_product_name,sum_by_gene_name)

#setup####
# filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
# attach(filtered_annotated)

#don't need no hypothetical proteins
#filtered_annotated <- filtered_annotated[grep(".*hypothetical protein.*",product_name,perl=T,invert=T),]
#do not actually use this
#sum_by_product_name <- read.csv("sum_by_product_name.csv")
#colnames(sum_by_product_name)[1]<-"product"
#attach(sum_by_product_name)
#
#don't need no hypothetical proteins
#sum_by_product_name <- sum_by_product_name[grep(".*hypothetical protein.*",product_name,perl=T,invert=T),]
#sum_by_product_name<-sum_by_product_name[,1:5]
#sum_by_product_name$product <- tolower(sum_by_product_name$product)
#sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
#colnames(sum_by_gene_name)[1]<-"gene"
#sum_by_gene_name<-sum_by_gene_name[,1:5]
#sum_by_gene_name$gene <- tolower(sum_by_gene_name$gene)
#
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

# patric_annotation <- read.delim("product_to_pathway.tab")
# patric_annotation$product <- tolower(patric_annotation$product)
# #get rid of quotes and %2c to be consistent
# gsub('\"','',patric_annotation$product)->patric_annotation$product
# gsub('%2c',',',patric_annotation$product)->patric_annotation$product
# print("making with_pathways")
# print(colnames(filtered_annotated))
# print(colnames(patric_annotation))
# with_pathways<-merge(x=filtered_annotated,y=patric_annotation,by.x="product_name",by.y="product")
# rm(patric_annotation,filtered_annotated)

####Start here again####

#get rid of blanks
# with_pathways<-with_pathways[grep(".+",with_pathways$pathway),]

#separate into individual pathways
# with_pathways<-separate_rows(with_pathways, pathway, sep = ";")
# with_pathways_no_dup<-with_pathways[!duplicated(with_pathways),]
# rm(with_pathways)

#and lets get our lps, polyamine and butyrate right now
# all_lps<-with_pathways_no_dup[grep(".*lipopolysaccharide biosynthesis.*",with_pathways_no_dup$pathway,perl = T,ignore.case = T),]
# write.table(all_lps,"all_lps_products.tab",sep = "\t", quote = T,row.names = F)
# rm(all_lps)

# all_butyrate<-with_pathways_no_dup[grep(".*Butanoate metabolism.*",with_pathways_no_dup$pathway,perl = T,ignore.case = T),]
# write.table(all_butyrate,"all_butyrate_products.tab",sep = "\t", quote = T,row.names = F)
# rm(all_butyrate)

# all_polyamine<-with_pathways_no_dup[grep(".*Arginine and proline metabolism.*",with_pathways_no_dup$pathway,perl = T,ignore.case = T),]
# write.table(all_polyamine,"all_polyamine_products.tab",sep = "\t", quote = T,row.names = F)
# rm(all_polyamine)

#add 'em up (this obviously causes over-estimation of true expression so it's all relative)
# sum_by_kegg_pathway<-rowsum(with_pathways_no_dup[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")],group = with_pathways_no_dup$pathway)
# sum_by_kegg_pathway<-sum_by_kegg_pathway[!is.na(sum_by_kegg_pathway$S1_FPM),]
# sum_by_kegg_pathway$sum<-rowSums(sum_by_kegg_pathway[,c("S1_FPM","S2_FPM","S3_FPM","S4_FPM")])
#
# #formatting for stupid bubble plot, which i bet could be done in R
# sum_by_kegg_pathway$Name<-row.names(sum_by_kegg_pathway)
#
# shortened<-sum_by_kegg_pathway[sum_by_kegg_pathway$sum>mean(sum_by_kegg_pathway$sum),]
# shortened<-shortened[order(shortened$sum , decreasing = T),]
# shortened<-shortened[,c("Name","S3_FPM","S4_FPM","S1_FPM","S2_FPM")]
# colnames(shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")
#
# write.table(shortened,"combined_sum_by_kegg_pathway_above_mean.tab", sep = "\t", quote = T,row.names = F)

#can do the rest interactively i think
