#save list of all rhizobiales MAGs
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)


# 1. load data

##mag info
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")

#rhizobialel MAGs
rhizo<-mags_role %>% filter(bac_order=="Rhizobiales") %>% 
  select(Genome) %>% distinct()
rhizo_fasta<-data.frame(paste0("analysis/05_MAGs/MAGs/bacs/",rhizo$Genome,".fa.gz"))
write.table(rhizo_fasta, "analysis/09_rhizobiales_phylogeny/list_rhizobiales_mags.txt", sep='\t',quote = F, row.names = F, col.names = F)
