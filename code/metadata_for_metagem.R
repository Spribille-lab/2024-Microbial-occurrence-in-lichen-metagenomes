#prepare metadata for metaGEM

## 1. misc
library(tidyverse)
library(stringr)
source("code/utils.R")

## 2. load data
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))
mags_role <-mags_role %>% left_join(checkm)  

### 3. Prepare metadata

selected<-mags_role %>% filter(metagenome=="X11" | metagenome=="VT1") %>%
  mutate(domain=ifelse(lineage_broad=="bacteria_other","bacteria","eukaryotes")) %>%
  mutate(phylum=ifelse(domain=="bacteria",bac_phylum,as.character(lineage_broad))) %>%
  select(Genome,metagenome,Lichen.metagenomes,depth_cov,compelteness,contamination,domain,phylum,bac_order,bac_family,bac_genus2,bat_bacteria)

colnames(selected)<-c("Genome","metagenome_ID","metagenomes_sample_name","depth_coverage","completness","contamination","domain","phylum","order","family","genus","full_taxonomy_bacteria")  
write.table(selected,"analysis/10_metaGEM/lichen_metadata_metaGEM.txt",sep="\t",quote = F, row.names = F)


## 4. prepared list of 15 chloro and 15 cyano lichens
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F)
colnames(depth_df)<-c("metagenome","depth")

###number of mags in all metagenomes
mag_number<-mags_role %>% group_by(metagenome, confirmed_role) %>% summarize(n_genome=n()) %>%
  pivot_wider(names_from=confirmed_role,values_from = n_genome, values_fill = 0)

###combine all data together
df<-depth_df %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  select(metagenome, Lichen.metagenomes, photobiont,depth) %>%
  left_join(mag_number)

###filter
df2<-df %>% filter(mycobiont==1,bacteria_other>0)

###trebouxioid lichens, ranged based on depth
treb<-df2 %>% filter(photobiont=="trebouxioid") %>% arrange(desc(depth)) %>% head(15)

###cyano lichens, ranged based on depth
cyan<-df2 %>% filter(photobiont=="cyano") %>% arrange(desc(depth)) %>% head(15)

### combine 
df3<-rbind(treb,cyan) %>% select(-mycobiont_missassigned,-cephalodia_cyano,-algae_other)
write.table(df3,"analysis/10_metaGEM/selected_mtg.txt",sep="\t",quote = F, row.names = F)

### list all mags in the selected metagenomes
df4<-mags_role %>% filter(metagenome %in% df3$metagenome) 
write.table(df4,"analysis/10_metaGEM/mags_in_selected_mtg.txt",sep="\t",quote = F, row.names = F)

         