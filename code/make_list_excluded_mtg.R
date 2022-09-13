##make a table with excluded metagenomes
#setwd("~/Documents/coverage")
library(tidyverse)
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

excl<-mags_role %>% filter(confirmed_role=="mycobiont_missassigned") %>%
  select(metagenome,Lichen.metagenomes,note)
write.table(excl,"analysis/05_MAGs/tables/excluded_metagenomes.txt",sep="\t",quote = F, row.names = F)
