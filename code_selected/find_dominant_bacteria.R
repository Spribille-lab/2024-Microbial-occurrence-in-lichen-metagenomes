#Identifying most common bacterial lineages

## 1. misc
library(tidyverse)
library(stringr)
source("code/utils.R")

## 2. load data
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mtg_info<-read.delim("results/all_metagenome_reanalysis.txt")

#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

#get taxonomy for bacterial mags 
bac_tax<-mags_role %>% filter(!is.na(mags_role$bac_phylum)) %>% select(Genome,bac_phylum,bac_order,bac_family,bac_genus,bac_genus2) %>% distinct()

###bwa read coverage
genome_coverage = read_csv("analysis/05_MAGs/tables/read_mapping/bwa_coverage.csv") 



## 3. Which bacteria occur most commonly, if we count multiple occurences of the same MAG?
### proccess coverage data
occurrences_genomes<-genome_coverage %>% pivot_longer(-Genome,names_to = "metagenome", values_to="coverage") %>%
  mutate(presence=ifelse(coverage>50,1,0)) %>% group_by(Genome) %>%
  summarise(occurrence=sum(presence))
### add taxonomy info 
occurrences_genomes<-occurrences_genomes %>% left_join(bac_tax) %>% drop_na()

### summarize by MAG
mag_occur<-occurrences_genomes %>% group_by(Genome) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
mag_occur<-mag_occur %>% left_join(bac_tax %>% select(Genome,bac_genus2,bac_family,bac_order,bac_phylum) %>%distinct())
write.table(genus_occur,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_mags_frequency.tsv",sep="\t",quote = F, row.names = F)

sum(mag_occur$occ_total[1:119])*100/sum(mag_occur$occ_total)
#119 top MAGs account for 50%

### summarize by genus
genus_occur<-occurrences_genomes %>% group_by(bac_genus2) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
genus_occur<-genus_occur %>% left_join(bac_tax %>% select(bac_genus2,bac_family,bac_order,bac_phylum) %>%distinct())
write.table(genus_occur,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_genera_frequency.tsv",sep="\t",quote = F, row.names = F)

sum(genus_occur$occ_total[1:12])*100/sum(genus_occur$occ_total)
#12 top genera account for 50%

### by family
family_occur<-occurrences_genomes %>% group_by(bac_family) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
family_occur<-family_occur %>% left_join(bac_tax %>% select(bac_family,bac_order,bac_phylum) %>%distinct())
write.table(family_occur,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_families_frequency.tsv",sep="\t",quote = F, row.names = F)

sum(family_occur$occ_total[1:3])*100/sum(family_occur$occ_total)
#3 top families account for 55%

###by order
order_occur<-occurrences_genomes %>% group_by(bac_order) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
order_occur<-order_occur %>% left_join(bac_tax %>% select(bac_order,bac_phylum) %>%distinct())
write.table(order_occur,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_order_frequency.tsv",sep="\t",quote = F, row.names = F)

sum(order_occur$occ_total[1:3])*100/sum(order_occur$occ_total)
#3 top orders account for 57%

## 4. Which bacteria have expansions in regrads with the number of diff MAGs?
### summarize by genus
genus_diversity<-bac_tax %>% group_by(bac_genus2) %>% summarise(n=n()) %>% arrange(desc(n))
genus_diversity<-genus_diversity %>% left_join(bac_tax %>% select(bac_genus2,bac_family,bac_order,bac_phylum) %>%distinct())
write.table(genus_diversity,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_genera_diversity.tsv",sep="\t",quote = F, row.names = F)

sum(genus_diversity$n[1:19])*100/sum(genus_diversity$n)
#19 top genera account for 50%

### by family
family_diversity<-bac_tax %>% group_by(bac_family) %>% summarise(n=n()) %>% arrange(desc(n))
family_diversity<-family_diversity %>% left_join(bac_tax %>% select(bac_family,bac_order,bac_phylum) %>%distinct())

write.table(family_diversity,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_families_diversity.tsv",sep="\t",quote = F, row.names = F)

sum(family_diversity$n[1:4])*100/sum(family_diversity$n)
#4 top families account for 50%

### by order
order_diversity<-bac_tax %>% group_by(bac_order) %>% summarise(n=n()) %>% arrange(desc(n))
order_diversity<-order_diversity %>% left_join(bac_tax %>% select(bac_order,bac_phylum) %>%distinct())
write.table(order_diversity,"analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_order_diversity.tsv",sep="\t",quote = F, row.names = F)

sum(order_diversity$n[1:4])*100/sum(order_diversity$n)
#4 top orders account for 53%

#in the dominant orders, what genera drive the abundance?
occurrences_genomes %>% filter(bac_order=="Rhizobiales") %>% group_by(bac_genus2) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
occurrences_genomes %>% filter(bac_order=="Acetobacterales") %>% group_by(bac_genus2) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
occurrences_genomes %>% filter(bac_order=="Acidobacteriales") %>% group_by(bac_genus2) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))
occurrences_genomes %>% filter(bac_order=="Sphingomonadales") %>% group_by(bac_genus2) %>% summarise(occ_total=sum(occurrence)) %>% arrange(desc(occ_total))


## 5. Are the top taxa the same between statistics?

###top genera
genus_occur$occ_rank<-1:nrow(genus_occur)
genus_diversity$diversity_rank<-1:nrow(genus_diversity)

genus_ranking<-genus_occur %>% left_join(genus_diversity) %>% select(bac_genus2, occ_rank,read_fraction_rank,diversity_rank,  bac_family, bac_order, bac_phylum)
##more similar

###top families
family_occur$occ_rank<-1:nrow(family_occur)
family_diversity$diversity_rank<-1:nrow(family_diversity)

family_ranking<-family_occur %>% left_join(family_diversity) %>% select(bac_family, occ_rank,read_fraction_rank,diversity_rank, bac_order, bac_phylum)
### for top 4, occurrense and reads are identical, diversity has same families in a different order

###top orders
order_occur$occ_rank<-1:nrow(order_occur)
order_diversity$diversity_rank<-1:nrow(order_diversity)

order_ranking<-order_occur %>% left_join(order_diversity) %>% select(bac_order, occ_rank,read_fraction_rank,diversity_rank,  bac_phylum)
### for top 4, same orders in a different order



## 6. make a list of dominant bacteria by photobiont type
#filter: remove metagenomes that didn't yeild mycobiont mag
mags_role_filtered<-mags_role %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont"]) %>%
  left_join(mtg_info,by=c("metagenome"="Run"))


trent_genus<-mags_role_filtered %>%
   filter(photobiont=="trentepohlioid", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

trent_family<-mags_role_filtered %>% 
  filter(photobiont=="trentepohlioid", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
           group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


treb_genus<-mags_role_filtered %>% 
  filter(photobiont=="trebouxioid", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

treb_family<-mags_role_filtered %>% 
  filter(photobiont=="trebouxioid", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

cyano_genus<-mags_role_filtered %>% 
  filter(photobiont=="cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

cyano_family<-mags_role_filtered %>% 
  filter(photobiont=="cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

## 11. make a list of dominant bacteria by mycobiont 

lichino_genus<-mags_role_filtered %>% 
  filter(class=="Lichinomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

lichino_family<-mags_role_filtered %>% 
  filter(class=="Lichinomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())



arthonio_genus<-mags_role_filtered %>% 
  filter(class=="Arthoniomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

arthonio_family<-mags_role_filtered %>% 
  filter(class=="Arthoniomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


dothideo_genus<-mags_role_filtered %>% 
  filter(class=="Dothideomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

dothideo_family<-mags_role_filtered %>% 
  filter(class=="Dothideomycetes", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


umb_genus<-mags_role_filtered %>% 
  filter(order=="Umbilicariales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

umb_family<-mags_role_filtered %>% 
  filter(order=="Umbilicariales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

calic_genus<-mags_role_filtered %>% 
  filter(order=="Caliciales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

calic_family<-mags_role_filtered %>% 
  filter(order=="Caliciales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

pert_genus<-mags_role_filtered %>% 
  filter(order=="Pertusariales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

pert_family<-mags_role_filtered %>% 
  filter(order=="Pertusariales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

baeom_genus<-mags_role_filtered %>% 
  filter(order=="Baeomycetales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

baeom_family<-mags_role_filtered %>% 
  filter(order=="Baeomycetales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

lecid_genus<-mags_role_filtered %>% 
  filter(order=="Lecideales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

lecid_family<-mags_role_filtered %>% 
  filter(order=="Lecideales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

telo_genus<-mags_role_filtered %>% 
  filter(order=="Teloschistales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

telo_family<-mags_role_filtered %>% 
  filter(order=="Teloschistales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


#verruc
ver_family<-mags_role_filtered %>% 
  filter(order=="Verrucariales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  group_by(bac_family) %>% summarise(occ_total=n()) %>% 
  arrange(desc(occ_total)) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


#peltigerales
peltigerales_cyano_genus<-mags_role_filtered %>% 
  filter(order=="Peltigerales",photobiont=="cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

#peltigerales chloro and cyano
peltigerales_cyano_genus<-mags_role_filtered %>% 
  filter(order=="Peltigerales",photobiont=="cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())

peltigerales_treb_genus<-mags_role_filtered %>% 
  filter(order=="Peltigerales",photobiont=="trebouxioid_cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct())


peltigerales_treb_family<-mags_role_filtered %>% 
  filter(order=="Peltigerales",photobiont=="trebouxioid_cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_family) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

peltigerales_cyano_family<-mags_role_filtered %>% 
  filter(order=="Peltigerales",photobiont=="cyano", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_family) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

mags_role_filtered %>% 
  filter(order=="Peltigerales", confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_family) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())

mags_role_filtered %>% 
  filter(photobiont %in% c("trebouxioid","trentepohlioid" ), confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_family) %>%
  summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>%
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct())


## 7. tie together previous section to create the supplementery tables
#### make a group as a single variable, a combination of taxonomy on diff. level and phtobiont
mags_role_filtered2<-mags_role_filtered %>%
  mutate(metagenome_group=ifelse(class=="Lecanoromycetes",order, class)) %>%
  mutate(metagenome_group=ifelse(order=="Peltigerales",paste0("Peltigerales_",photobiont), metagenome_group)) %>% 
  filter(metagenome_group!="Unknown") 

explore_groups<-mags_role_filtered2 %>% group_by(metagenome_group) %>% summarize(n=n())

#### list groups that have bacteria in them
groups_with_bacteria<-mags_role_filtered2 %>% filter(confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
  select(metagenome_group) %>% distinct() 

#### make table of dominant genera by group
top_bac_genus<-function(group_of_int){
  df<- mags_role_filtered2 %>% 
    filter(metagenome_group==group_of_int, confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% group_by(bac_genus2) %>%
    summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>% head(10) %>%
    mutate(metagenome_group=group_of_int)
  return(df)
}

l<-lapply(groups_with_bacteria$metagenome_group %>% sort,top_bac_genus)
genus_join_table<-do.call(rbind,l) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus) %>% distinct()) %>%
  select(metagenome_group, bac_genus2,bac_family,bac_order,bac_phylum,occ_total)

#### make table of dominant families by group
top_bac_fam<-function(group_of_int){
  df<- mags_role_filtered2 %>% 
    filter(metagenome_group==group_of_int, confirmed_role %in% c("bacteria_other","cephalodia_cyano", "photobiont_cyano")) %>% 
    group_by(bac_family) %>%
    summarise(occ_total=n()) %>% arrange(desc(occ_total)) %>% head(10) %>%
    mutate(metagenome_group=group_of_int)
  return(df)
}

l2<-lapply(groups_with_bacteria$metagenome_group %>% sort,top_bac_fam)
fam_join_table<-do.call(rbind,l2) %>% 
  left_join(bac_tax %>% select(-Genome,-bac_genus,-bac_genus2) %>% distinct()) %>%
  select(metagenome_group, bac_family,bac_order,bac_phylum,occ_total)

#### add dominant across all dataset
genus_occur_top<-genus_occur %>% mutate(metagenome_group="Whole dataset") %>%
  head(10) %>% select(metagenome_group,bac_genus2,bac_family,bac_order,bac_phylum,occ_total)
genus_join_table<- rbind(genus_occur_top, genus_join_table)

family_occur_top<-family_occur %>% mutate(metagenome_group="Whole dataset") %>%
  head(10) %>% select(metagenome_group,bac_family,bac_order,bac_phylum,occ_total)
fam_join_table<- rbind(family_occur_top, fam_join_table)

write.table(genus_join_table,"results/bacterial_genera_by_lichen_group.txt",sep="\t",quote = F, row.names = F)
write.table(fam_join_table,"results/bacterial_families_by_lichen_group.txt",sep="\t",quote = F, row.names = F)




## 8. Select high-quality MAGs from the 12 top-occurrence genera
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))
mag_occur <-mag_occur %>% left_join(checkm)  

good_mags<-mag_occur %>% filter(bac_genus2 %in% genus_occur$bac_genus2[1:12]) %>% 
  filter(compelteness>95,contamination<10) 
###make locus tags
first<-str_replace_all(good_mags$bac_genus2, "[^[:alnum:]]", "") %>% substr(1,8) %>% toupper()
second<-gsub('[a-zA-Z]+_(.*)_.*_.*', '\\1', good_mags$Genome) 
locustag<-paste0(first,second)

### added two more MAGs for annotation: the two top-occurrence Lichinhabitans MAGs 
mags_to_annotate<-good_mags$Genome


###added the only good VCDI01 (=Lichenicoccus) MAG
v<-mag_occur %>% filter(bac_genus2=="VCDI01",compelteness>95,contamination<10)
v$Genome

mags_to_annotate<-c(mags_to_annotate,v$Genome)
mags_to_annotate_files<-paste("analysis/05_MAGs/MAGs/bacs/",mags_to_annotate,".fa.gz",sep="")
write.table(data.frame(mags_to_annotate),"analysis/07_annotate_MAGs/mag_list.txt",sep="\t",quote = F, row.names = F,col.names = F)


## 9.  save table with info on selected mags
mag_info<-mag_occur  %>% select(-c(occ_total,occ_rang,strain_heterog,taxonomy)) %>% filter(Genome %in% mags_to_annotate)
write.table(mag_info,"analysis/07_annotate_MAGs/mag_table_info.txt",sep="\t",quote = F, row.names = F)




