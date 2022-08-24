##make a table with excluded metagenomes
#setwd("~/Documents/coverage")
library(tidyverse)
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")

excl<-mags_role %>% filter(confirmed_role=="mycobiont_missassigned") %>%
  select(metagenome,Lichen.metagenomes,note)
write.table(excl,"results/tables/excluded_metagenomes.txt",sep="\t",quote = F, row.names = F)
