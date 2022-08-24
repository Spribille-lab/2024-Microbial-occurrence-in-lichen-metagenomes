###check the presence of kegg modules of interest in the 65 selected MAGs

setwd("~/Documents/coverage")
library(tidyverse)
source("code/utils.R")

#load data
module_df<-read.table("analysis/07_annotate_MAGs/summarized_outputs/kegg_module_matrix.txt",header  = T)

mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_family2=ifelse(bac_family=="Unknown",paste(bac_order," fam.",sep=""),bac_family))
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family2," gen. sp.",sep=""),bac_genus))
mag_taxonomy <- mags_role %>% select(Genome,bac_genus2,bac_family2,bac_order,bac_phylum) %>% distinct()

annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
annotated_mags<-annotated_mags %>% left_join(mag_taxonomy)

##remove the two lower completness Lichenihabitans MAGs
annotated_mags2<- annotated_mags %>%
  filter(Genome!="private_T1916_metawrap_bin.6" & Genome != "private_T1889_metawrap_bin.7")


##load table with module description
module_descr<-read.table("analysis/07_annotate_MAGs/summarized_outputs/kegg_module_description.txt",header  = T,sep="\t")

modules_of_interest<-module_descr %>% 
  filter(grepl("Cobalamin",Definition)  | grepl("Biotin",Definition)  | grepl("Riboflavin",Definition) | grepl("Thiamine",Definition) )

#join tables
df<-module_df %>% pivot_longer(-Genome,names_to="Module",values_to="presence") %>%
  filter(Module %in% modules_of_interest$Module & Genome %in% annotated_mags2$Genome) %>%
  left_join(module_descr) %>% left_join(mag_taxonomy)
df2<-df %>% group_by(bac_genus2,Module) %>% summarize(perc_pres=sum(presence)/n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=perc_pres) %>%
  left_join(module_descr)
write.table( df2, "analysis/07_annotate_MAGs/summarized_outputs/modules_of_interest_freq.txt", sep='\t',quote = F, row.names = F, col.names = T)

