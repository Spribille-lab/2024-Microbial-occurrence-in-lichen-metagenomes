##misc: a script that calculates some of the stats used in the text
source("code/utils.R")
library(tidyverse)
#percentage of mycobiont mags among all fungal mags
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")


mags_role %>% filter(lineage_broad == "Fungi") %>% group_by(confirmed_role) %>%
  summarize(n=n()) %>% pivot_wider(names_from = confirmed_role,values_from = n) %>%
  mutate(myco_percent = (mycobiont + mycobiont_missassigned)/(mycobiont + mycobiont_missassigned + cypho + fungi_other + tremella))

#what do we know about the taxonomy of fungi_other based on the BAT assignments?
mags_role %>% filter(confirmed_role == "fungi_other") %>% group_by(BAT_euk) %>%
  summarize(n=n())

#are there any metagenomes that don't have a paul-prooduced assembly?
assembly_list<-read.delim("manuscript/tmp_files/n",header=F)
mtg_list<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
mtg_list %>% filter(!(Run %in% assembly_list$V1))

#what portion of bacterial mags comes from proteobacteria and acidobacteria
bacteria<-read.delim("analysis/05_MAGs/trees/bac_itol/itol_gtdb-layer.txt",sep=",",skip = 4,header=F)
bac_summ<-bacteria %>% group_by(V3) %>% summarize(n=n()) %>% pivot_wider(names_from = V3,values_from = n) %>%
  mutate(proteo_acido=(Proteobacteria+Acidobacteriota)/nrow(bacteria))
bac_summ$proteo_acido

#how many bacterial mags were detected in more than one metagenome?
mags_role %>% filter(lineage_broad == "bacteria_other" | lineage_broad == "Cyanobacteria" ) %>% group_by(Genome) %>%
  summarize(n=n()) %>%arrange(desc(n)) %>% mutate(binary=ifelse(n>1,T,F)) %>% 
  group_by(binary) %>% summarize(n=n())

mags_role %>% filter(lineage_broad == "bacteria_other" | lineage_broad == "Cyanobacteria" ) %>% group_by(Genome) %>%
  summarize(n=n()) %>%arrange(desc(n)) 

 

#metagenomes per bacterial mag
mtg_count<-mags_role %>% filter(lineage_broad == "bacteria_other" | lineage_broad == "Cyanobacteria" ) %>% 
  group_by(Genome) %>% summarize(n=n())  %>% arrange(desc(n)) 
mtg_count %>% mutate(more5=ifelse(n<=5,"less5","more5")) %>% group_by(more5) %>% summarise(n=n())
ggplot(mtg_count, aes(x=n))+geom_histogram()

  #separate by group
mtg_count<-mags_role %>% filter(lineage_broad == "bacteria_other" | lineage_broad == "Cyanobacteria" ) %>% 
  group_by(Genome,bac_family) %>% summarize(n=n())  %>% arrange(desc(n)) 

fam_to_plot<-c("Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","Beijerinckiaceae")
ggplot(mtg_count %>% filter(bac_family %in% fam_to_plot), aes(x=n))+
  geom_histogram()+
  facet_wrap(.~bac_family)

#algal mags number of ocurrences
mags_role %>% filter(confirmed_role== "photobiont_chloro") %>% group_by(Genome) %>%
  summarize(n=n())  %>% arrange(desc(n)) 

#what are the metagenomes with lichenihabitans?
mags_role %>% filter(bac_genus=="Lichenihabitans") %>% dplyr::select(metagenome,Lichen.metagenomes) %>%
  unique() 

mags_role %>% filter(bac_genus=="Lichenihabitans") %>% group_by(Genome) %>%
  summarise(n=n())

#how many times the two dominant mags co-occurred
mags_role %>% filter(Genome=="private_T1889_metawrap_bin.7"| Genome=="public_SRR14722130_metawrap_bin.2") %>%
  select(metagenome,Genome,depth_cov)%>%
  group_by(Genome) %>% pivot_wider(names_from=Genome, values_from=depth_cov,values_fill=0) %>%
  filter(private_T1889_metawrap_bin.7>0 & public_SRR14722130_metawrap_bin.2>0)




#lichenihabitans mag occurrences per study
lichenihab_source<-mags_role %>% left_join(mtg_list,by=c("metagenome"="Run")) %>% 
  filter(bac_genus=="Lichenihabitans") %>%
  group_by(Genome,Source) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = Source, values_from = n,values_fill=0) %>%
  mutate()

acido_source<-mags_role %>% left_join(mtg_list,by=c("metagenome"="Run")) %>% 
  filter(bac_order=="Acidobacteriales") %>%
  group_by(Genome,Source) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = Source, values_from = n,values_fill=0) %>%
  mutate()

aceto_source<-mags_role %>% left_join(mtg_list,by=c("metagenome"="Run")) %>% 
  filter(bac_order=="Acetobacterales") %>%
  group_by(Genome,Source) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = Source, values_from = n,values_fill=0) %>%
  mutate()

#what metagenomes had the yeasts
bp<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F,col.names=c("metagenome","bp"))

mags_role %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
   left_join(bp) %>% 
   filter(confirmed_role=="tremella" | confirmed_role=="cypho") %>%
   select(Genome, metagenome, Lichen.metagenomes.y,confirmed_role,depth_cov,bp,order,Source)                       

#select photobionts from letharias to see if there are different in each metagenome
mags_role %>% left_join(mtg_info,by=c("metagenome"="Run")) %>%
  filter(Source=="Tuovinen et al. 2019" & confirmed_role == "photobiont_chloro") %>%
  select(Genome, metagenome,Lichen.metagenomes.y,confirmed_role)

#how many metagenomes had an algal mag based on the bwa alignments?
mags_role %>% filter(confirmed_role=="photobiont_chloro") %>% 
  select(metagenome) %>% unique() %>% nrow()


#what fraction of all bacterial MAG occurences come from the top 4 families?
fam_to_plot<-c("Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","Beijerinckiaceae")
(mags_role %>% filter(bac_family %in% fam_to_plot) %>% nrow())/ (mags_role %>% filter(lineage_broad %in% c("bacteria_other","Cyanobacteria")) %>% nrow())

#what fraction of all MAGs come from the top 4 families?
fam_to_plot<-c("Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","Beijerinckiaceae")
mags_in_fam<-mags_role %>% filter(bac_family %in% fam_to_plot) %>% select(Genome) %>% distinct()
all_bac_mags<-mags_role %>% filter(lineage_broad %in% c("bacteria_other","Cyanobacteria")) %>% select(Genome) %>% distinct()
nrow(mags_in_fam)/nrow(all_bac_mags)


