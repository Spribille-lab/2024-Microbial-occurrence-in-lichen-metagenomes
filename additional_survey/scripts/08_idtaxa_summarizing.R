# Script for summarizing the idtaxa and metaxa2 outputs for bacterial families

# adapted from Gulnara Tagridzhanova's script

# load in necessary packages
source("scripts/utils.R")
library(tidyverse)

##make a table with mag occurrences of the groups of interest
groups_to_include<-c("Trebouxiophyceae",
                     "Acidobacteriaceae",
                     "Beijerinckiaceae",
                     "Acetobacteraceae",
                     "Sphingomonadaceae",
                     "Cystobasidiomycetes",
                     "Tremellomycetes")



# 2. get info from metaxa for 18S of CYpho, Tremella, and Trebouxia. screening of assemblies and reads
matrix_level5<-read.delim("results/metaxa_level_5_combined.txt")

matrix_level5_long<-matrix_level5 %>% 
  pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% 
  summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% 
  mutate(occurrence=ifelse(abundance>0,1,0)) %>% 
  filter(!str_detect(metagenome, 'SRR10277352|SRR19097886'))

##summarize by group, reduce the table only to target groups 
metaxa_occ_target<-matrix_level5_long %>% 
  dplyr::filter(occurrence>0, 
                grepl( "Cystobasidiomycetes", Taxa) |
                  grepl( "Tremellomycetes", Taxa) |
                  grepl( "Trebouxi", Taxa)) %>%
  mutate(label = case_when(str_detect(Taxa, 'Cystobasidiomycetes') ~
                             'Cystobasidiomycetes',
                           str_detect(Taxa, 'Beijerinckiaceae') ~ 
                             'Beijerinckiaceae',
                           str_detect(Taxa, 'Tremellomycetes') ~ 
                             'Tremellomycetes',
                           str_detect(Taxa, 'Trebouxi') ~ 
                             'Trebouxiophyceae',
                           TRUE ~ NA_character_)) %>%
  dplyr::select(-Taxa,-abundance) %>% 
  group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% 
  mutate(occurrence=ifelse(n>0,1,0))

# 3. get info from idtaxa for 16S screening
idtaxa <- read.delim("results/idtaxa/idtaxa_reads_all.txt") %>% 
  mutate(metagenome = str_remove(file, 
                                  'results/metaxa_bacteria/') %>% 
           str_remove('.bacteria.fasta'),
         type = 'reads') %>% 
  separate(col = assignment, 
           into = c('drop_root',
                    'drop_kingdom',
                    'drop_kingdom2',
                    'bac_phylum',
                    'drop_class',
                    'bac_order',
                    'label',
                    'drop_genus'),
           sep = ';',
           fill = 'right',
           remove = FALSE) %>% 
  dplyr::select(-contains('drop_')) %>% 
  filter(!str_detect(metagenome, 'SRR10277352|SRR19097886'))

# these two should not be included
# SRR10277352 # taken from culture
# SRR19097886 # Sarea


# metagenomes_to_keep  # ADD HERE

# idtaxa <- idtaxa %>% filter(metagenome %in% metagenomes_to_keep) #added filtering step to make sure that only intended metagenomes are counted, i.e. the metagenomes with a mycobiont mag and those that were not excluded as misidentified

#get family_level assignment


idtaxa_occ_target<-idtaxa %>%
  dplyr::filter(label %in% groups_to_include) %>%
  group_by(type, metagenome, label) %>%
  summarize(n=n()) %>% 
  mutate(occurrence=ifelse(n>0,1,0))

metaxa_occ_target %>% pull(metagenome) %>% unique() %>% length()
idtaxa_occ_target %>% pull(metagenome) %>% unique() %>% length()

# 4. get together metaxa and idtaxa
rdna <- rbind(idtaxa_occ_target,metaxa_occ_target)  %>%
  pivot_wider(id_cols = -n,
              values_from = occurrence, 
              values_fill = 0, 
              names_from = label)


# write.table(rdna,"analysis/tables/occurrence_rDNA_groups_of_interest.tsv",sep="\t",quote = F, row.names = F)


# 8. calculate stats
n_mtg <- rdna %>% pull(metagenome) %>% unique() %>% length()
n_mtg


##make a list of top bacterial families in assemblies and reads. ranked by the % of metagenomes they occurred in

#reads
top_bacteria_reads<-idtaxa %>% 
  dplyr::filter(type=="reads",!is.na(label)) %>% 
  dplyr::select(metagenome,label) %>%
  dplyr::distinct() %>%  
  dplyr::group_by(label) %>%
  dplyr::summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(percentage=n*100/n_mtg) %>% # ADD HERE
  dplyr::arrange(dplyr::desc(n))


# write.table(top_bacteria_reads,"analysis/tables/idtaxa_reads_bac_families.tsv",sep="\t",quote = F, row.names = F)

#combine both tables to make a supplementary table
top_idtaxa<- top_bacteria_reads %>% 
  mutate(Source="Reads") %>%
  mutate(Number_of_metagenomes_detected = n,
                 Prevalence = n*100/n_mtg) %>% 
  left_join(idtaxa %>% 
              select(label,bac_order,bac_phylum) %>%
              distinct()) %>%
          select(Source,label,bac_order,bac_phylum,
                 Prevalence,Number_of_metagenomes_detected) %>% 
  mutate(label = str_replace(label, '_', ' ') %>% 
           str_to_title()) %>% 
  dplyr::select(-Source) %>% 
  head(25) %>% 
  write.csv(file = 'output/top_25_bacteria_publication.csv', 
            row.names = FALSE,
            na = '')

# write.table(top_idtaxa,"analysis/tables/idtaxa_top_bac_families.tsv",sep="\t",quote = F, row.names = F)


# Lichenihabitans

metaxa_occ_target

idtaxa %>% 
  separate(col = assignment, 
           into = c('root',
                    'kingdom',
                    'kingdom2',
                    'phylum',
                    'class',
                    'order',
                    'family',
                    'genus'),
           sep = ';',
           fill = 'right',
           remove = FALSE) %>% 
  dplyr::select(metagenome, genus) %>% 
  group_by(metagenome) %>% 
  summarize(Lichenihabitans = if_else(any(str_detect(genus, 'Lichenihabitans')), 
                                      'yes',
                                      'no') %>% 
              replace_na('no')) %>%
  group_by(Lichenihabitans) %>% 
  count()

print('Percent of metagenomes that contain Lichenihabitans:')
164/(79+164)



# Analysis for all groups of interest
prev_focal <- rdna %>% 
  pivot_longer(cols = Acetobacteraceae:Tremellomycetes) %>% 
  group_by(name) %>% 
  summarize(n = sum(value),
            percent = n/rdna %>% pull(metagenome) %>% unique() %>% length()) %>% 
  arrange(desc(percent))

prev_focal

write.csv(prev_focal, 'output/new_record_percentages_publication.csv', row.names = FALSE, na = '')


