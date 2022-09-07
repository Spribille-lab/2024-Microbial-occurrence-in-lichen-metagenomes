#Creating combined kegg mapper files by genus

## 1. misc
library(tidyverse)
library(stringr)
library(vegan)
source("code/utils.R")

## 2. read data
mags_role<-read.delim(results/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")

## 3. get taxonomy
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())

###how many mags per genus?
table(annotated_mags$bac_genus2)

## 4. read kegg files
read_kegg<-function(mag_list){

read_kegg_by_mag<-function(mag){
  filename<-paste0("analysis/07_annotate_MAGs/",mag,"/",mag,".kegg.mapper.txt")
  keggfile<-read.delim(filename,header=F)
  name_df<-data.frame(V1=paste0("# ",mag),V2="")
  keggfile<-rbind(name_df,keggfile)
  return(keggfile)}

l<-lapply(mag_list,read_kegg_by_mag)

return(l)   
}

## 5. aplied the function to genera 
read_kegg_genus<-function(genus){
mags_from_genus<-annotated_mags$Genome[annotated_mags$bac_genus2==genus]
list<-read_kegg(mags_from_genus)
filename<-paste0('analysis/07_annotate_MAGs/summarized_outputs/',genus,'.kegg.combined.txt')
lapply(list, function(x) write.table( data.frame(x), filename, append= T, sep='\t',quote = F, row.names = F, col.names = F))
}

for (genus in unique(annotated_mags$bac_genus2)){
read_kegg_genus(genus)
}

## 6. Created a table for all MAGs
read_kegg_table<-function(mag_list){
  
  read_kegg_by_mag<-function(mag){
    filename<-paste0("analysis/07_annotate_MAGs/",mag,"/",mag,".kegg.mapper.txt")
    keggfile<-read.delim(filename,header=F,na.strings = "")
    keggfile$Genome<-mag
    return(keggfile)}
  
  l<-lapply(mag_list,read_kegg_by_mag)
  
  return(l)   
}

l<-read_kegg_table(annotated_mags$Genome)
kegg_combined<-do.call(rbind,l)
write.table( kegg_combined, "analysis/07_annotate_MAGs/summarized_outputs/selected_mags_kegg_combined.txt", sep='\t',quote = F, row.names = F, col.names = T)






