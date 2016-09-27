setwd("/Users/Scott/Google Drive/Hurwitz Lab/cuffdiff_reallyall")

filtered_annotated <- read.csv("diff_exp_for_all_bact.csv")
sum_by_product_name <- read.csv("sum_by_product_name.csv")
sum_by_gene_name <- read.csv("sum_by_gene_name.csv")

mschaedleri <- read.csv("MschaedleriFeatures.csv")
mschaedleri_exp <- merge(mschaedleri,filtered_annotated,by.x="RefSeq.Locus.Tag",by.y="tracking_id",all=F)
mschaedleri_exp_vs_all_products <- merge(mschaedleri_exp,sum_by_product_name,by.x="product_name",by.y="X",all=F)

write.csv(mschaedleri_exp_vs_all_products,"mschaedleri_contribution_to_set.csv")
