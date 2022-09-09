#setwd("~/Documents/coverage")
library(tidyverse)
library(taxize)
library(myTAI)

#get final list of MAGs
mags_full_table<-read.csv("analysis/05_MAGs/tables/mag_sizes.csv",header=F)
colnames(mags_full_table)<-c("fasta","size")

#get coverage
cov<-read.csv("analysis/05_MAGs/tables/cmseq_coverage.csv")
cov$fasta<-gsub('cov_', '', cov[,1])
cov<-cov[,2:3]

#get taxonomy estimates
##from bacteria
bat_bacteria<-read.delim("analysis/05_MAGs/tables/gtdbtk.bac120.summary.tsv")[,1:2]
colnames(bat_bacteria)<-c("mag","bat_bacteria")

## from euks
busco<-read.csv("analysis/05_MAGs/tables/eukcc_busco_table.csv")[,c(1,10)]

bat_euk<-read.delim("analysis/05_MAGs/tables/BAT_euks.bin2classification.txt")[,c(1,4)]
bat_euk$taxid<-gsub('.*;', '', bat_euk[,2])
###transform taxids into names
taxon_summary<-ncbi_get_taxon_summary(bat_euk$taxid)
bat_euk$BAT_euk<-taxon_summary[,2]
bat_euk<-bat_euk %>% select(X..bin,lineage,BAT_euk)
colnames(bat_euk)<-c("fasta","lineage","BAT_euk")

#join tables

mag_taxonomy<-left_join(mags_full_table,busco) %>% left_join(bat_euk)  %>% left_join(cov) 
mag_taxonomy$mag<-gsub('.{3}$', '', mag_taxonomy$fasta) #remove '.fa' from the MAG names
mag_taxonomy<-left_join(mag_taxonomy,bat_bacteria) 

write.table(mag_taxonomy,"analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv",sep="\t",quote = F, row.names = F)


#prepare table to use as supplementary
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","completeness","contamination","strain_heterog","taxonomy"))
mags_full_table$Genome<-str_replace(mags_full_table$fasta,".fa","")
bac<-checkm %>% inner_join(bat_bacteria,by=c("Genome"="mag")) %>% 
  left_join(mags_full_table) %>% mutate(classification=bat_bacteria) %>%
  select(Genome,completeness,contamination,size,classification) 

eukcc<-read.csv("analysis/05_MAGs/tables/eukcc_busco_table.csv",header=T)
euk<-eukcc %>% left_join(size) %>% mutate(classification=BUSCO_lineage) %>%
  select(Genome,completeness,contamination,size,classification)

write.table(rbind(bac,euk),"analysis/05_MAGs/tables/MAG_suppl_table.tsv",sep="\t",quote = F, row.names = F)



