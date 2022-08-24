#checking psudofinder results
##are loci identified as pseudogenes located on contig ends?

#BiocManager::install("GenomicRanges")
#BiocManager::install("GenomicFeatures")

library(ape)
library(Biostrings)
source("code/utils.R")
library(tidyverse)

##mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% dplyr:::select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())


#read an example of pseudofinder gff
gff<-read.gff("analysis/07_annotate_MAGs/pseudofinder/private_T1916_metawrap_bin.6_pseudofinder_pseudos.gff", na.strings = c(".", "?"), GFF3 = TRUE)

assembly<-readDNAStringSet("analysis/07_annotate_MAGs/private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.fna")
names(assembly)<-str_replace(names(assembly),"gnl\\|UoN\\|","")
assembly_info<-data.frame("seqid"=names(assembly),"length"=width(assembly))

gff<-gff %>% left_join(assembly_info) %>% mutate(dist_to_end=length-end)
nrow(gff)

gff_ends<-gff %>% filter(start < 10 | dist_to_end < 10 )
nrow(gff_ends)

#make a function to check all MAGs annotated with pseudofinder
pseudo_ann_mags<-list.files("analysis/07_annotate_MAGs/pseudofinder/", "*_pseudofinder_pseudos.gff", all.files =F , full.names =T)

get_marginal_counts<-function(gff_path){
  gff<-read.gff(gff_path, na.strings = c(".", "?"), GFF3 = TRUE)
  mag_name<- gff_path %>% str_replace("analysis/07_annotate_MAGs/pseudofinder//","") %>%
    str_replace("_pseudofinder_pseudos.gff","")
  assembly_path<-paste0("analysis/07_annotate_MAGs/",mag_name,"/",mag_name,".fna")
  assembly<-readDNAStringSet( assembly_path)
  names(assembly)<-str_replace(names(assembly),"gnl\\|UoN\\|","")
  assembly_info<-data.frame("seqid"=names(assembly),"length"=width(assembly))
  gff<-gff %>% left_join(assembly_info) %>% mutate(dist_to_end=length-end)
  total<-nrow(gff)
  gff_ends<-gff %>% filter(start < 10 | dist_to_end < 10 )
  marginal<-nrow(gff_ends)
  return(marginal/total)
  }
  
l<-lapply(pseudo_ann_mags,get_marginal_counts)
out<-do.call(rbind,l)
max(out)

#get info from the pseudofinder log files
pseudo_logs<-list.files("analysis/07_annotate_MAGs/pseudofinder/", "*_pseudofinder_log.txt", all.files =F , full.names =T)
logfile<-pseudo_logs[1]

get_pseudo_table<-function(logfile){
  mag_name<- logfile %>% str_replace("analysis/07_annotate_MAGs/pseudofinder//","") %>%
    str_replace("_pseudofinder_log.txt","")
  f<-read.delim(logfile,skip=30,header=F)
  table<-f[c(1:3,5:11),]
  table$V1<-str_replace_all(table$V1, "[[:punct:]]", "")
  table$V1<-str_replace_all(table$V1, " ", "_")
  table$V2<-as.character(table$V2) %>% as.numeric()
  table<-table %>% pivot_wider(names_from = V1, values_from=V2)
  table$Genome<-mag_name
  return(table)}

l<-lapply(pseudo_logs,get_pseudo_table)
df<-do.call(rbind,l)
df <- df %>% mutate(pseudogenes_perc=Pseudogenes_total*100/(Pseudogenes_total+Intact_genes))
df <- annotated_mags %>% left_join(df)

df %>% group_by(bac_genus2) %>% 
  summarize(avg_pseudo=mean(Pseudogenes_total,na.rm=T),avg_pseudo_perc=mean(pseudogenes_perc,na.rm=T),
            avg_gene_original=mean(Initial_ORFs,na.rm=T)) 

#get a list of pseudogenes
get_psudo_lists<-function(gff_path){
  gff<-read.gff(gff_path, na.strings = c(".", "?"), GFF3 = TRUE)
  mag_name<- gff_path %>% str_replace("analysis/07_annotate_MAGs/pseudofinder//","") %>%
    str_replace("_pseudofinder_pseudos.gff","")
  gff$old_locus_tag<-gsub("^.*old_locus_tag=(.*)$","\\1",gff$attributes)
  locus_list<-data.frame("locustag"=gff$old_locus_tag)
  locus_list$attributes<-gff$attributes
  locus_list<-locus_list%>% filter(!grepl("_ign_", locustag)) %>%
    separate_rows(1,sep = ",")
  locus_list$Genome<-mag_name
  locus_list$type <- "pseudogene"
  return(locus_list)
}

l<-lapply(pseudo_ann_mags,get_psudo_lists)
pseudo_list<-do.call(rbind,l)
write.table(pseudo_list,"analysis/07_annotate_MAGs/summarized_outputs/pseudogene_list.txt",sep="\t",quote = F, row.names = F)


## what are KOs of the pseudogenes?
# are any of genes from KO of interest pseudogenized?
kegg_combined<-read.delim("analysis/07_annotate_MAGs/summarized_outputs/kegg_combined.txt")
kegg_of_interest<-read.delim("analysis/07_annotate_MAGs/kegg_of_interest.txt")
kegg_list<-kegg_of_interest[,1]


pseudo_df<-kegg_combined %>% filter(V2 %in% kegg_list) %>% filter(V1 %in% pseudo_list$locustag) %>%
  left_join(kegg_of_interest,by=c("V2"="KO"))


#are those pseudogenes the only members of their family in a given MAG? or are there any intact genes?
pseudo_df2<-pseudo_df %>% mutate(genome_KO=paste(Genome,V2,sep="_"))


kegg_combined2 <- kegg_combined %>% left_join(pseudo_list,by=c("V1"="locustag","Genome"="Genome")) %>%
  mutate(count=ifelse(is.na(type),1,0)) %>% mutate(genome_KO=paste(Genome,V2,sep="_"))

lost<-kegg_combined2 %>% filter(genome_KO %in% pseudo_df2$genome_KO) %>% 
  group_by(Genome,V2) %>% summarize(count=sum(count)) %>%
  filter(count==0) %>% left_join(kegg_of_interest,by=c("V2"="KO")) %>%
  left_join(annotated_mags)
write.table(lost,"analysis/07_annotate_MAGs/summarized_outputs/pseudogenized_kegg_of_interest.txt",sep="\t",quote = F, row.names = F)



# 5. What are the KOs of all pseudogenized genes? what are the KOs that are lost to pseudogenes?
pseudo_df_full<-kegg_combined  %>% left_join(pseudo_list,by=c("V1"="locustag","Genome"="Genome")) %>%
  filter(!is.na(V2)) %>% mutate(count=ifelse(is.na(type),1,0)) %>%
 mutate(genome_KO=paste(Genome,V2,sep="_"))

pseudo_df_sum<-pseudo_df_full %>% filter(genome_KO %in% pseudo_df_full$genome_KO[pseudo_df_full$count==0]) %>% 
  group_by(Genome,V2) %>% summarize(count=sum(count)) %>%
  left_join(annotated_mags) %>%
  filter(count==0) %>%
  group_by(bac_family,V2) %>%
  summarize(n=n()) %>%
left_join(kegg_of_interest,by=c("V2"="KO"))

#looking at all pseudogenes with kegg assignments
pseudo_kegg_ann<-pseudo_df_full %>% filter(count==0) %>% left_join(kegg_of_interest,by=c("V2"="KO"))


 