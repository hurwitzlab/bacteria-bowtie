setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffnorm-out")

filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
sum_by_product_name <- read.csv("sum_by_product_name.csv")
colnames(sum_by_product_name)[1]<-"product"
sum_by_gene_name <- read.csv("sum_by_gene_name.csv")
colnames(sum_by_gene_name)[1]<-"gene"

#for LPS pathway####
lps_path<-read.table("LPS_search_list",sep = "\t",comment.char = "#")
lowercase_lps<-data.frame(tolower(lps_path[,1]))
lps_path<-lowercase_lps
colnames(lps_path)<-"V1"

just_lps_products<-merge(x=sum_by_product_name,by.x="product",y=lps_path,by.y="V1",all=F)
just_lps_genes<-merge(x=sum_by_gene_name,by.x="gene",y=lps_path,by.y="V1",all=F)

polyamines<-read.table("polyamine-search-list",sep="\t",comment.char = "#")
polyamines<-data.frame(tolower(polyamines[,1]))
colnames(polyamines)<-"V1"
polyamines<-unique(polyamines)

just_poly_products<-merge(x=sum_by_product_name,by.x="product",y=polyamines,by.y="V1",all=F)
just_poly_genes<-merge(x=sum_by_gene_name,by.x="gene",y=polyamines,by.y="V1",all=F)

#for m schaedleri####
mschaedleri <- read.csv("MschaedleriFeatures.csv")
mschaedleri_exp <- merge(mschaedleri,filtered_annotated,by.x="RefSeq.Locus.Tag",by.y="tracking_id",all=F)
mschaedleri_exp_vs_all_products <- merge(mschaedleri_exp,sum_by_product_name,by.x="product_name",by.y="X",all=F)

write.csv(mschaedleri_exp_vs_all_products,"mschaedleri_contribution_to_set.csv")
