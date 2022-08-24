#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(ape)
library(phytools)

# read bacterial tree
bac_tree<-read.newick("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")
old_labels<-data.frame(bac_tree$tip.label)



#make new labels
bat_bacteria<-read.delim("analysis/05_MAGs/tables/gtdbtk.bac120.summary.tsv")[,1:2]
bat_bacteria$classification<-gsub("[][!#$%()*,.:;<=>@^_`|~.{} ]", "", bat_bacteria$classification)
bat_bacteria$new_name<-paste(bat_bacteria$user_genome,bat_bacteria$classification,sep='_')

lables_both<-left_join(old_labels,bat_bacteria,by=c("bac_tree.tip.label"="user_genome"))
bac_tree$tip.label<-lables_both$new_name
write.tree(bac_tree, file='analysis/05_MAGs/trees/renamed/gtdbtk.bac120.user_msa.fasta_renamed.treefile')


#make new iTOL annotation file
itol_old<-read.csv("analysis/05_MAGs/trees/bac_itol/itol_gtdb-layer.txt",skip=4,header=F)

itol_new<-left_join(itol_old,lables_both,by=c("V1"="bac_tree.tip.label"))

itol_new2<-itol_new[c(5,2,3)]
write.table(itol_new2,"analysis/05_MAGs/trees/renamed/itol_gtdb-layer_renamed.txt",sep=",",quote = F, row.names = F, col.names=F)


## make iTOL annotation for the number of occurrences per mag

mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

occurrences<-mags_role %>% filter(breadth>=50) %>%
  group_by(Genome) %>% summarize(occurrences=n()) 

occurrences_df<-itol_new %>% left_join(occurrences,by=c("V1"="Genome")) %>% select(new_name,occurrences)


cat("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Occurrences\nCOLOR,#16537e\nDATA\n",file="analysis/05_MAGs/trees/renamed/itol_gtdb-layer_renamed_occurrences.txt")
write.table(occurrences_df,"analysis/05_MAGs/trees/renamed/itol_gtdb-layer_renamed_occurrences.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)




