#setwd("~/Documents/gulya/coverage")
library(tidyverse)

#read table for SRA runs from PRJNA700635 (Genome streamlining)
sra_table1<-read.csv('analysis/03_lendemer_reanalysis/SraRunInfo_PRJNA700635.csv')

#read table for SRA runs from PRJNA731936 (Lichens of the Southern Appalachians)
sra_table2<-read.csv('analysis/03_lendemer_reanalysis/SraRunInfo_PRJNA731936.csv')

#combine
sra_table<-rbind(sra_table1,sra_table2)

#extracte FEN number
sra_table$FEN_id<-sub("([^_]+_+[^_]+).*", "\\1",sra_table$SampleName )



#read rable from Lendemer et al.
table<-read.csv('analysis/03_lendemer_reanalysis/ajb21339-sup-0002-appendixs2.csv',sep='\t')
table$FEN_id<-sub("([^_]+_+[^_]+).*", "\\1",table$FEN )

#join tables saving only enries that exist in both tables
combined_table<-inner_join(sra_table,table) %>% 
  select(Run,SampleName,ScientificName,Current_Determination,Basidiomycete.Hit.) %>%
  arrange(Run) %>% mutate(processed=ifelse(Run %in% sra_table1$Run,T,F))

#manually check all entries where the species name is not identical between the SRA table and table from Lendemer
combined_table$ScientificName<-as.character(combined_table$ScientificName)
combined_table$Current_Determination<-as.character(combined_table$Current_Determination)
combined_table %>% filter(ScientificName!=Current_Determination)

#save SRA ids
write.table(combined_table$Run,"analysis/03_lendemer_reanalysis/sra_ids_lendemer.txt",col.names = F,row.names = F, quote = F)

#save the table
write.table(combined_table,"analysis/03_lendemer_reanalysis/lendemer_table_sra.txt",col.names = T,row.names = F, quote = F,sep="\t")

