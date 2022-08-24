#setwd("~/Documents/coverage")
library(tidyverse)
library(ape)
library(phytools)

# read the tree
#euk_tree<-read.newick("analysis/05_MAGs/trees/euk_tree_with_ref/lichen_mags_iqtree.tree")
#euk_tree<-read.newick("analysis/05_MAGs/trees/euk_tree_with_ref/lichen_mags.contree")
#euk_tree<-read.newick("analysis/05_MAGs/trees/euk_tree_with_ref/lichen_mags.treefile")
euk_tree<-read.newick("analysis/05_MAGs/trees/eukaryotes.treefile")

old_labels<-data.frame(euk_tree$tip.label)


#get info on the MAGs from this study
mag_taxonomy<-read.delim("analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv")
#mag_taxonomy<-read.delim("analysis/05_MAGs/tables/MAG_taxonomy_and_role.tsv")
mag_taxonomy$Run<-gsub('[a-zA-Z]+_(.*)_.*_.*', '\\1', mag_taxonomy$mag)
metagenome_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
mag_taxonomy<-mag_taxonomy %>% left_join(metagenome_info) 
mag_taxonomy$new_name<-paste(mag_taxonomy$fasta,mag_taxonomy$BAT_euk,mag_taxonomy$Lichen.metagenomes,sep="_")

#get the info on the reference genomes
ref<-read.csv("analysis/05_MAGs/trees/reference_genomes/reference_genomes_full_table.csv")
ref$new_name<-paste(ref$ID,ref$Species,sep="_")

#get all info into one dataset
both_labels_mags<-mags_taxonomy %>% select(fasta, new_name) %>% rename(old_name=fasta)
both_labels_ref<-ref %>% select(ID,new_name)%>% rename(old_name=ID)
both_labels<-rbind(both_labels_mags,both_labels_ref)

both_labels<-left_join(old_labels,both_labels,by=c("euk_tree.tip.label"="old_name"))

#annotate the tree
euk_tree$tip.label<-both_labels$new_name
write.tree(euk_tree, file='analysis/05_MAGs/trees/eukaryotes_renamed.treefile')


#iterative renaming - after checking the tree, readjust the labels
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

mags_role2<-mags_role %>% select(Genome, confirmed_role) %>% 
  group_by(Genome, confirmed_role) %>%summarize(n=n()) %>% 
  pivot_wider(values_from=n,names_from=confirmed_role,values_fill=0)
mags_role2<-mags_role %>% select(Genome, metagenome,confirmed_role) %>%
  left_join(metagenome_info %>% select(Run,Lichen.metagenomes),by=c("metagenome"="Run"))


#if a mag labeled as a mycobiont, put it in it's name all lichens it's assigned as mycbiont of
mycobiont_list<-mags_role2 %>% filter(confirmed_role=="mycobiont" | confirmed_role=="mycobiont_missassigned") %>%
  group_by(Genome) %>% summarize(all_mtg=paste(Lichen.metagenomes, collapse=" "))

mycobiont_list$all_mtg<-paste("mycobiont",mycobiont_list$all_mtg,sep=" ")

#if a fungal MAG is classified as fungi_other and IS NOT a mycobiont in another metagenome
fungi_other_list<-mags_role2 %>% filter(confirmed_role=="fungi_other" & !(Genome %in% mycobiont_list$Genome)) %>%
  group_by(Genome) %>% summarize(all_mtg=paste(Lichen.metagenomes, collapse=" "))
fungi_other_list$all_mtg<-paste("fungi_other",fungi_other_list$all_mtg,sep=" ")


fungal_labels<-rbind(mycobiont_list,fungi_other_list)
fungal_labels$iterative_label<-paste(fungal_labels$Genome, fungal_labels$all_mtg)
fungal_labels$Genome<-paste(fungal_labels$Genome,".fa",sep='')

iterative_label<-both_labels %>% left_join(fungal_labels,by=c("euk_tree.tip.label"="Genome"))

#for the others (ref. genomes and algae) keep the same labels
iterative_label$iterative_label[is.na(iterative_label$iterative_label)]<-as.character(iterative_label$new_name[is.na(iterative_label$iterative_label)])
iterative_label$iterative_label[iterative_label$euk_tree.tip.label=="public_SRR12240178_metabat2_bin.6.fa"]<-"public_SRR12240178_metabat2_bin.6 mycobiont Xanthoparmelia spp"


#rename the tree
euk_tree2<-euk_tree
euk_tree2$tip.label<-iterative_label$iterative_label
write.tree(euk_tree2, file='analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile')



