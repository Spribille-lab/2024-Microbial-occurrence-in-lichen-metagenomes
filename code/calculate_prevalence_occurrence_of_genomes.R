#Count occorrences of different types of MAGs

## 1. misc
library(tidyverse)
library(stringr)
source("code/utils.R")

## 2. load data
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")

#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

## 3. For key groups, count % of MAG occurrences they represent
total_occ<-mags_role %>% nrow()
total_bac_occ_number<-mags_role %>% filter(lineage_broad %in% c("bacteria_other", "Cyanobacteria")) %>% nrow

total_bac_occ_number/total_occ #% of bacterial occurrence among all occurrences

acetobacteraceae_occ<- mags_role %>% filter(bac_family == "Acetobacteraceae") %>% nrow
acidobacteriaceae_occ<- mags_role %>% filter(bac_family == "Acidobacteriaceae") %>% nrow
beijerinckiaceae_occ<- mags_role %>% filter(bac_family == "Beijerinckiaceae") %>% nrow
sphingomonadaceae_occ<- mags_role %>% filter(bac_family == "Sphingomonadaceae") %>% nrow
###% of occurrences of the 4 most common bacterial families
(acetobacteraceae_occ+acidobacteriaceae_occ+beijerinckiaceae_occ+sphingomonadaceae_occ)/total_occ #of all occurrences
(acetobacteraceae_occ+acidobacteriaceae_occ+beijerinckiaceae_occ+sphingomonadaceae_occ)/total_bac_occ_number #of bacterial occurrences

fungal_occ<-mags_role %>% filter(lineage_broad == "Fungi") %>% nrow
fungal_occ/total_occ #% of fungal occurrence among all occurrences

algal_occ<-mags_role %>% filter(lineage_broad == "Chlorophyta") %>% nrow
algal_occ/total_occ #% of chlorophyte occurrence among all occurrences

cyano_occ<-mags_role %>% filter(lineage_broad == "Cyanobacteria") %>% nrow
cyano_occ/total_occ #% of cyanobacterial occurrence among all occurrences

## 4. For key groups, count % of metagenomes they occurred in
### make a list of metagenomes that yeilded at least one MAG
metagenomes_with_mags<-mags_role %>% select(metagenome) %>% distinct()

### % of metagenomes that contained at least one MAG from the 4 families
four_families_mtg<- mags_role %>% 
  filter(bac_family %in% c("Acetobacteraceae", "Acidobacteriaceae", "Beijerinckiaceae", "Sphingomonadaceae")) %>%
  select(metagenome) %>% distinct() %>% nrow()
four_families_mtg/nrow(metagenomes_with_mags)

### % of metagenomes that contained at least one fungal MAG 
fungi_mtg<- mags_role %>% 
  filter(lineage_broad == "Fungi") %>%
  select(metagenome) %>% distinct() %>% nrow()
fungi_mtg/nrow(metagenomes_with_mags)

### % of metagenomes that contained at least one algal MAG 
alga_mtg<- mags_role %>% 
  filter(lineage_broad == "Chlorophyta") %>%
  select(metagenome) %>% distinct() %>% nrow()
alga_mtg/nrow(metagenomes_with_mags)

### % of metagenomes that contained at least one Cyano MAG 
cyano_mtg<- mags_role %>% 
  filter(lineage_broad == "Cyanobacteria") %>%
  select(metagenome) %>% distinct() %>% nrow()
cyano_mtg/nrow(metagenomes_with_mags)



