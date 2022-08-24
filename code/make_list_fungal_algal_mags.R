#setwd("~/Documents/coverage")
#made a list of fungal and algal mags for the phylogenomic trees

library(tidyverse)
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

fungal_mag<-mags_role %>% filter(lineage_broad=="Fungi") %>% 
  mutate(filename_full=paste(fasta,".gz",sep="")) %>% select(filename_full) %>%
  unique() 
write.table(fungal_mag,"analysis/05_MAGs/tables/list_fungal_MAGs.txt",sep="\t",quote = F, col.names= F, row.names = F)

algal_mag<-mags_role %>% filter(lineage_broad=="Chlorophyta") %>% 
  mutate(filename_full=paste(fasta,".gz",sep="")) %>% select(filename_full) %>%
  unique() 
write.table(algal_mag,"analysis/05_MAGs/tables/list_algal_MAGs.txt",sep="\t",quote = F, col.names= F, row.names = F)

#make a list of reference genomes by taxonomy
ref<-read.delim("analysis/05_MAGs/trees/reference_genomes/reference_genomes_full_table.csv")
algal_ref<-ref %>% filter(grepl("Chlorophyta",Taxonomy))
fungal_ref<-ref %>% filter(grepl("mycota",Taxonomy))
