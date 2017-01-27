#To be used after "pathway_mapping_on_HPC.R"

setwd("~/bacteria-bowtie/scripts/R-interactive/")

shortened<-read.table("combined_sum_by_kegg_pathway_above_mean.tab",header = T)
more_shortened<-shortened[1:30,]
colnames(more_shortened)=c("Name","S+H- (Control)","S-H- (SMAD3 Knockout)","S+H+ (H. hepaticus only)","S-H+ (Combined)")

write.table(more_shortened,"combined_sum_by_kegg_pathway_above_mean.tab", sep = "\t", quote = T,row.names = F)

#Have to run this in an external shell cuz ... perl... grumble
#system("./bubble.sh combined_sum_by_kegg_pathway_above_mean.tab CombinedBubble")

system("cp CombinedBubble.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")

more_shortened$combined_effect<-log(more_shortened$`S-H+ (Combined)`/more_shortened$`S+H- (Control)`)
more_shortened$Hhep_effect<-log(more_shortened$`S+H+ (H. hepaticus only)`/more_shortened$`S+H- (Control)`)
more_shortened$smad_effect<-log(more_shortened$`S-H- (SMAD3 Knockout)`/more_shortened$`S+H- (Control)`)
more_shortened<-more_shortened[order(more_shortened$combined_effect , decreasing = T),]

effects <- data.matrix(more_shortened[,6:8])

row.names(effects)<-more_shortened$Name

x=effects

oldPar <- par(no.readonly = T)

myColors=colorRampPalette(c("Blue","Yellow"))

#a heatmap, cuz why not!
heatmap(x, Rowv=NA, Colv=NA, col = myColors(255),scale="none", margins=c(5,5), cexCol=1, labCol = c("Combined", "H. hepaticus","SMAD3-KO"))

new_bubble_source <- more_shortened[,1:5]

write.table(new_bubble_source,"sum_by_kegg_pathway_ordered_by_combined_effect.tab", sep = "\t", quote = T,row.names = F)

system("source ~/.bash_profile && ./bubble.sh sum_by_kegg_pathway_ordered_by_combined_effect.tab CombinedBubble_orderedbyeffect")

system("cp CombinedBubble_orderedbyeffect.pdf '/Users/Scott/Google Drive/Hurwitz Lab/manuscripts/'")
