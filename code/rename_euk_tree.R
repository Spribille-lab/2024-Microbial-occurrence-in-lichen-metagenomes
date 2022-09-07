#setwd("~/Documents/coverage")
library(tidyverse)
library(ape)
library(phytools)

# fungal tree

## 1. read the tree
fun_tree<-read.newick("analysis/05_MAGs/trees/eukaryotes/Fungi/phylogeny/Concatenated_IQTREE/concat.contree")
fun_old_labels<-data.frame(fun_tree$tip.label)
colnames(fun_old_labels)<-"old_name"
metagenome_info<-read.delim(results/tables/all_metagenome_reanalysis.txt")
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_putative_roles_bwa.tsv")


## 2. get the info on the reference genomes
ref<-read.csv("results/tables/reference_genomes_full_table.csv")
ref$new_name<-paste(ref$ID,ref$Species,sep="_")

## 3. connect names in the ref genome table to the names in the tree tips

### get genome IDs out of the tip names
ref_genome_ids<-fun_old_labels %>% filter(grepl("GC",old_name))
ref_genome_ids$ID<-sub("^([^_]*_[^_]*).*", "\\1", ref_genome_ids$old_name)

### add manually the ref genomes not from NCBI
rhomi<-c("Rhomi1_AssemblyScaffolds","Rhomi1")
cypho<-c("Cyphobasidialed_GT57","CAJHEP01")
trem<-c("Tremellales_GT57","GCA_904859935.1")
ref_genome_ids<-rbind(ref_genome_ids,rhomi,cypho,trem)

## 4. make new names for reference genomes
ref_genome_ids<-ref_genome_ids %>% left_join(ref) %>% mutate(new_name=paste(ID,Species,sep="_")) %>%
  select(old_name,new_name)

## 5. make new names for the MAGs by adding info
mags_role2<-mags_role %>% select(Genome, metagenome,putative_role) %>%
  left_join(metagenome_info %>% select(Run,Lichen.metagenomes),by=c("metagenome"="Run"))

##for each MAG make a list of metagenomes where it is a putative mycobiont
mycobiont_list<-mags_role2 %>% filter(putative_role=="mycobiont") %>%
  group_by(Genome) %>% summarize(all_mtg=paste(Lichen.metagenomes, collapse=" "))

mycobiont_list$all_mtg<-paste("mycobiont",mycobiont_list$all_mtg,sep=" ")

#if a fungal MAG is classified as fungi_other and IS NOT a mycobiont in another metagenome
fungi_other_list<-mags_role2 %>% filter(putative_role=="fungi_other" & !(Genome %in% mycobiont_list$Genome)) %>%
  group_by(Genome) %>% summarize(all_mtg=paste(Lichen.metagenomes, collapse=" "))
fungi_other_list$all_mtg<-paste("fungi_other",fungi_other_list$all_mtg,sep=" ")

##prep data for renaming
fungal_mags_labels<-rbind(mycobiont_list,fungi_other_list) %>% 
  mutate(new_name=paste0(Genome,all_mtg,sep="_"),old_name=Genome) %>%
  select(old_name,new_name)
###manually changed one mag name, since it was way too long
fungal_mags_labels$new_name[fungal_mags_labels$old_name=="public_SRR12240178_metabat2_bin.6"]<-"public_SRR12240178_metabat2_bin.6 mycobiont Xanthoparmelia spp"



## 6. rename tree tips

both_labels<-rbind(fungal_mags_labels,ref_genome_ids)
both_labels<-left_join(fun_old_labels,both_labels)

fun_tree$tip.label<-both_labels$new_name
write.tree(fun_tree, file='analysis/05_MAGs/trees/eukaryotes/Fungi/phylogeny/Concatenated_IQTREE/concat_putative_renamed.contree')


# algal tree: here I only renamed reference genomes, and kept MAG names as is
## 1. read the tree
alg_tree<-read.newick("analysis/05_MAGs/trees/eukaryotes/Algae/phylogeny/Concatenated_IQTREE/concat.contree")
alg_old_labels<-data.frame(alg_tree$tip.label)
colnames(alg_old_labels)<-"old_name"

## 2. connect names in the ref genome table to the names in the tree tips

### get genome IDs out of the tip names
ref_genome_ids<-alg_old_labels %>% filter(grepl("GC",old_name))
ref_genome_ids$ID<-paste0(sub("^([^_]*_[^_]*).*", "\\1", ref_genome_ids$old_name),".1")

## 3. add new info
ref_genome_ids<-ref_genome_ids %>% left_join(ref) %>% mutate(new_name=paste(ID,Species,Taxonomy,sep="_")) %>%
  select(old_name,new_name)

## 4. rename tree tips

both_labels<-left_join(alg_old_labels,ref_genome_ids)
both_labels$new_name[is.na(both_labels$new_name)]<-both_labels$old_name[is.na(both_labels$new_name)]

alg_tree$tip.label<-both_labels$new_name
write.tree(alg_tree, file='analysis/05_MAGs/trees/eukaryotes/Algae/phylogeny/Concatenated_IQTREE/concat_renamed.contree')


