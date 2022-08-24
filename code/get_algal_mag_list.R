#get the list of good algal mags

setwd("~/Documents/coverage")


## misc
library(tidyverse)

##load info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
qc<-read.delim("analysis/05_MAGs/tables/eukcc_busco_table.csv",sep=",")

##combine info and select only "good" algal MAGs  based on their EukCC scores
table<-mags_role %>% left_join(qc) %>% filter(BUSCO_lineage=="chlorophyta_odb10") %>%
  select(Genome,BAT_euk,completeness,contamination) %>% distinct() %>%
  filter(completeness > 90,contamination<10)
  
##save the list
data.frame(table$Genome)
write.table(data.frame(table$Genome), "analysis/11_algal_MAGs/good_algal_mags_list.txt", sep='\t',quote = F, row.names = F, col.names = F)

