library(tidyverse)

## Select high-quality MAGs from the 13 top-occurrence genera
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))
mag_occur<-read.delim("analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_mags_frequency.tsv")
mag_occur <-mag_occur %>% left_join(checkm)  

good_mags<-mag_occur %>% filter(bac_genus2 %in% genus_occur$bac_genus2[1:13]) %>% 
  filter(compelteness>95,contamination<10) 
mags_to_annotate<-good_mags$Genome
write.table(data.frame(mags_to_annotate),"analysis/07_annotate_MAGs/mag_list.txt",sep="\t",quote = F, row.names = F,col.names = F)

mags_to_annotate_files<-paste("analysis/05_MAGs/MAGs/bacs/",mags_to_annotate,".fa.gz",sep="")
mags_to_annotate_df<-data.frame(mags_to_annotate,mags_to_annotate_files)
colnames(mags_to_annotate_df)<-c("mag","mag_file_path")
write.table(mags_to_annotate_df,"analysis/07_annotate_MAGs/mag_table.txt",sep="\t",quote = F, row.names = F)

##save table with info
mag_info<-mag_occur  %>% select(-c(occ_total,strain_heterog,taxonomy)) %>% filter(Genome %in% mags_to_annotate)
write.table(mag_info,"results/tables/selected_mag_info.txt",sep="\t",quote = F, row.names = F)

