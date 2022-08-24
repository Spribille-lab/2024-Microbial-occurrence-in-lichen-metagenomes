#using multivar stats to find bacterial groups

setwd("~/Documents/coverage")


## 1. misc
library(tidyverse)
library(patchwork)
library(stringr)
library(vegan)
library(reshape2)
library(apcluster)
library(qualpalr)
library(DescTools) 
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(dbscan)
library(simplifyEnrichment)
library(ape)
library(phangorn)
require("phytools")
source("code/utils.R")
source("code/kegg_module_reconstruct.R")

## 2. read data
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

## 3. get taxonomy
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

mag_taxonomy <- mags_role %>% select(Genome,bac_genus2,bac_family,bac_order,bac_phylum) %>% distinct()

###list of bacterial mags
bac_mags<-mags_role[!(is.na(mags_role$bac_genus)),1] %>% unique()

###mag occurrence
mag_occur<-mags_role %>% group_by(Genome) %>% summarize(n=n())


## 4. Created a table for all MAGs
read_kegg_table<-function(mag_list){
  
  read_kegg_by_mag<-function(mag){
    filename<-paste0("analysis/07_annotate_MAGs/annotation_ALL_mags/",mag,".kegg.mapper.txt")
    keggfile<-read.delim(filename,header=F,na.strings = "")
    keggfile$Genome<-mag
    return(keggfile)}
  
  l<-lapply(mag_list,read_kegg_by_mag)
  
  return(l)   
}


l<-read_kegg_table(bac_mags)
kegg_combined<-do.call(rbind,l)
write.table( kegg_combined, "analysis/07_annotate_MAGs/summarized_outputs/all_mags_kegg_combined.txt", sep='\t',quote = F, row.names = F, col.names = T)
kegg_combined<-read.delim( "analysis/07_annotate_MAGs/summarized_outputs/all_mags_kegg_combined.txt")

## 5. Prepare data for nmds (excluded private_T1889_metawrap_bin.7 MAG, because it has low completeness)
matrix<-kegg_combined %>% select(Genome,V2) %>% 
  group_by(Genome,V2) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from=V2,values_from=n,values_fill=0) %>% data.frame()
rownames(matrix)<-matrix$Genome

nmds=metaMDS(matrix %>% select(-Genome),k=2,trymax=200)
#No convergence -- monoMDS stopping criteria:

## 6. Exploring ordinations

#Only retain MAGs that occur > 2 times
mag_occur2<-mag_occur %>% filter(n>2)

matrix2<-kegg_combined %>% filter(Genome %in% mag_occur2$Genome) %>% select(Genome,V2) %>% 
  group_by(Genome,V2) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from=V2,values_from=n,values_fill=0) %>% data.frame()
rownames(matrix2)<-matrix2$Genome

nmds=metaMDS(matrix2 %>% select(-Genome),k=2,trymax=200)

###plot

data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Genome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mag_taxonomy)  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$KO <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

order<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_order),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

phylum<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_phylum),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

order+phylum

fam<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Acetobacteraceae" = "red", "Acidobacteriaceae"="blue","Beijerinckiaceae"="brown","Sphimgomonadaceae"="violet","Nostocaceae"="#4CC5C2","UBA10450"="orange")) +
  coord_equal() +
  theme_bw()

## only retain MAGs that are >90% complete
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))
complete_mags<-checkm %>% filter(compelteness>90,contamination<10)

matrix3<-kegg_combined %>% filter(Genome %in% complete_mags$Genome) %>% select(Genome,V2) %>% 
  group_by(Genome,V2) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from=V2,values_from=n,values_fill=0) %>% data.frame()
rownames(matrix3)<-matrix3$Genome

nmds=metaMDS(matrix3 %>% select(-Genome),k=2,trymax=200)
#No convergence

## only retain MAGs that are >90% complete AND occurr > 2
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))
complete_mags<-checkm %>% filter(compelteness>90,contamination<10)

matrix4<-kegg_combined %>% filter(Genome %in% complete_mags$Genome & Genome %in% mag_occur2$Genome) %>% select(Genome,V2) %>% 
  group_by(Genome,V2) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from=V2,values_from=n,values_fill=0) %>% data.frame()
rownames(matrix4)<-matrix4$Genome

nmds=metaMDS(matrix4 %>% select(-Genome),k=2,trymax=200)

###plot

data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Genome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mag_taxonomy)  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$KO <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

order<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_order),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

phylum<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_phylum),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

order+phylum

fam<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Acetobacteraceae" = "red", "Acidobacteriaceae"="blue","Beijerinckiaceae"="brown","Sphimgomonadaceae"="violet","Nostocaceae"="#4CC5C2","UBA10450"="orange")) +
  coord_equal() +
  theme_bw()

## 7. Optimal ordination: first tried, with mags with >2 occurrences

mag_occur2<-mag_occur %>% filter(n>2)

matrix2<-kegg_combined %>% filter(Genome %in% mag_occur2$Genome) %>% select(Genome,V2) %>% 
  group_by(Genome,V2) %>% summarise(n=n()) %>% ungroup() %>%
  pivot_wider(names_from=V2,values_from=n,values_fill=0) %>% data.frame()
rownames(matrix2)<-matrix2$Genome

nmds=metaMDS(matrix2 %>% select(-Genome),k=2,trymax=200)

###plot

data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Genome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mag_taxonomy)  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$KO <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

order<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_order),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

phylum<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_phylum),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()



fam<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Acetobacteraceae" = "red", "Acidobacteriaceae"="blue","Beijerinckiaceae"="brown","Sphimgomonadaceae"="violet","Nostocaceae"="#4CC5C2","UBA10450"="orange")) +
  coord_equal() +
  theme_bw()

order+phylum+fam
ggsave("analysis/07_annotate_MAGs/summarized_outputs/nmds_all_mags_kegg.png",device="png",height = 6,width=20,bg="white")


## 8. Reconstruct kegg modules and use them for nmds. only used genomes that are >90% complete
km_dia<-KMdiagram_fetcher(create_RData=T, path="analysis/07_annotate_MAGs/rdata", new_date=Sys.Date())
km_matrix<-matrix3 %>% select(-Genome)
kmreco<-KMreco(indata=km_matrix, km_dia, len_breaks=NULL, allowed_gaps=c(1)) 
module_matrix<-kmreco[[4]]


### remove mags with low occurrence
module_df<-data.frame(module_matrix)
module_df$Genome<- matrix3$Genome
write.table( module_df, "analysis/07_annotate_MAGs/summarized_outputs/kegg_module_matrix.txt", sep='\t',quote = F, row.names = F, col.names = T)
module_df<-read.table("analysis/07_annotate_MAGs/summarized_outputs/kegg_module_matrix.txt",header  = T)

rownames(module_df)<-module_df$Genome
nmds_module=metaMDS(module_df %>% filter(Genome %in% mag_occur2$Genome) %>% select(-Genome), k=2,trymax=200)

###plot

data.scores <- as.data.frame(scores(nmds_module))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Genome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mag_taxonomy)  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds_module, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$KO <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

order<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_order),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

phylum<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_phylum),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

fam<-ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Acetobacteraceae" = "red", "Acidobacteriaceae"="blue","Beijerinckiaceae"="brown","Sphingomonadaceae"="violet","Nostocaceae"="#4CC5C2","UBA10450"="orange")) +
  coord_equal() +
  theme_bw()

order+phylum+fam
ggsave("analysis/07_annotate_MAGs/summarized_outputs/nmds_all_mags_kegg_module.png",device="png",height = 6,width=20,bg="white")


## 9. Hierarchical clustering - based on kegg modules
gnm_cor_module = cor(t(module_df %>% select(-Genome))) #correlation matrix

### 9a. selecting best clustering method

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison.pdf",width=10,height=7)
compare_clustering_methods(gnm_cor_module, plot_type = "heatmap")
dev.off()

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(gnm_cor_module,plot_type = "heatmap")
dev.off()
#### results: kmeans, hdbscan and apcluster look optimal

### 9b. clustering
set.seed(123)
hdbscan_gnm <- cluster_by_hdbscan(gnm_cor_module)
hdbscan_results<-data.frame(cbind(module_df$Genome,hdbscan_gnm))
colnames(hdbscan_results)<-c("Genome","hdbscan")

set.seed(123)
apcluster_gnm <- cluster_by_apcluster(gnm_cor_module)
apcluster_results<-data.frame(cbind(module_df$Genome,apcluster_gnm))
colnames(apcluster_results)<-c("Genome","apcluster")

set.seed(123)
kmeans_gnm <- cluster_by_kmeans(gnm_cor_module)
kmeans_results<-data.frame(cbind(module_df$Genome,kmeans_gnm))
colnames(kmeans_results)<-c("Genome","kmeans")


module_df2<-module_df %>% left_join(mag_taxonomy) %>% 
  select(Genome,bac_family,bac_order,bac_phylum) %>%
  left_join(hdbscan_results) %>%
  left_join(kmeans_results) %>%
  left_join(apcluster_results)
write.table(module_df2, "analysis/07_annotate_MAGs/summarized_outputs/kegg_modules_clustering.txt", sep='\t',quote = F, row.names = F, col.names = T)

### heatmap
ra = rowAnnotation(df=data.frame(bac_family = module_df2$bac_family,
                                 bac_order = module_df2$bac_order,
                                 bac_phylum = module_df2$bac_phylum,
                                 apcluster=module_df2$apcluster,
                                 kmeans=module_df2$kmeans,
                                 hdbscan=module_df2$hdbscan)
                   #col = list(type = type_color),
                   #annotation_legend_param = list(type=list(labels=type_ordered,at=type_ordered,title = "function"))
)

###this heatmap clusters rows and columns by hclust!
Heatmap(module_df %>% select(-Genome), show_column_names = F, show_row_names = F, name = " ",
        right_annotation = ra)

###same heatmap but with phylogenomic tree as a way to order rows
tip_to_remove<-bactree$tip.label[!(bactree$tip.label %in% module_df$Genome)]
bactree_filtered<-drop.tip(bactree, tip_to_remove)
write.tree(bactree_filtered, "analysis/05_MAGs/tmptree2.tre")
row_clust = ReadDendrogram("analysis/05_MAGs/tmptree2.tre")

pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_results.pdf",width=12,height=15)
Heatmap(module_df %>% select(-Genome), show_column_names = F, show_row_names = F, name = " ",
        #cluster_rows = row_clust,
        column_title = "Fig. S7. Functional clustering of bacterial MAGs. Each row represents a MAG, each column is a \nKEGG module. The bars on the left show taxonomic position of the MAG on three levels (family, \norder, and phylum) and the cluster it was assigned to by three clustering methods \n(apcluster, kmeans, and hdbscan). The dendrograms show the hierarchical clustering done using \nthe hclust function in R. For this analysis, we only used MAGs with at least 90% completeness.",
        right_annotation = ra)
dev.off()


### 9c  taxonomic coherence
### e.g. clustering vs taxonomy. also based on Zoccarato et al. 2022
###TC = Ngfc/Ntaxon. 
  #Ngfc = # of genomes included in the cluster. 
  #Ntaxon = # of genomes desendent from the last common ancestor of the cluster

bactree <- read.tree("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")

###function that takes a vector with node names and a tree and returns # of nodes descending from the MRCA
get_Ntaxon<-function(tip_list,bactree){
  mrca<-findMRCA(bactree, tip_list, type="node")
  l<-Descendants(bactree, mrca, type = "tips")[[1]]
  return(bactree$tip.label[l])
}

###function that process a cluster
process_cluster<-function(cluster_name,clustering_df,bactree){
  clustering_df$include<-clustering_df[,2]==cluster_name
  cluster<-clustering_df %>% filter(include==T)
  Ngfc<-cluster %>% nrow()
  if(Ngfc==1){
    Ntaxon<-1
  }else{
  Ntaxon<-get_Ntaxon(cluster$Genome,bactree) %>% data.frame() %>% 
    filter(. %in% clustering_df$Genome) %>% data.frame() %>% nrow() #only count mags that are in the clustering dataset (otherwise there are a lot of mags that were excluded from this analysis on early stages because of incompletenes)
  }
  out<-data.frame("cluster"=cluster_name,"Ngfc"=Ngfc,"Ntaxon"=Ntaxon)
  return(out)
}

l<-lapply(apcluster_results$apcluster %>% unique,process_cluster,clustering_df=apcluster_results,bactree=bactree)
aplcuster_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(hdbscan_results$hdbscan %>% unique,process_cluster,clustering_df=hdbscan_results,bactree=bactree)
hdbscan_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(kmeans_results$kmeans %>% unique,process_cluster,clustering_df=kmeans_results,bactree=bactree)
kmeans_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

## result: in all, the majority of clusters is not monophyletic. 
###to account for possibility of singletones, allowed exclusion of singletones
singletones_fam<-mag_taxonomy%>% filter(Genome %in% bac_mags) %>% group_by(bac_family) %>% summarize(n=n()) %>% filter(n==1)
singletones_mag<-mag_taxonomy %>% filter(bac_family %in% singletones_fam$bac_family)

process_cluster2<-function(cluster_name,clustering_df,bactree){
  clustering_df$include<-clustering_df[,2]==cluster_name
  cluster<-clustering_df %>% filter(include==T) %>% filter(!(Genome %in% singletones_mag$Genome))
  Ngfc<-cluster %>% nrow()
  if(Ngfc<2){
    Ntaxon<-1
  }else{
  Ntaxon<-get_Ntaxon(cluster$Genome,bactree) %>% data.frame() %>% 
    filter(. %in% clustering_df$Genome) %>% nrow() #only count mags that are in the clustering dataset (otherwise there are a lot of mags that were excluded from this analysis on early stages because of incompletenes)
  }
  out<-data.frame("cluster"=cluster_name,"Ngfc"=Ngfc,"Ntaxon"=Ntaxon)
  return(out)
}

l<-lapply(apcluster_results$apcluster %>% unique,process_cluster2,clustering_df=apcluster_results,bactree=bactree)
aplcuster_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(hdbscan_results$hdbscan %>% unique,process_cluster2,clustering_df=hdbscan_results,bactree=bactree)
hdbscan_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(kmeans_results$kmeans %>% unique,process_cluster2,clustering_df=kmeans_results,bactree=bactree)
kmeans_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)
##didn't improve much

### visualize: how namy taxonomic rank per cluster?
ntaxa_graph<-function(clustering_df,mag_taxonomy){
  clustering_df$cluster<-clustering_df[,2]
  clustering_df <- clustering_df %>% left_join(mag_taxonomy) %>% filter(cluster>0)

  n_order<-clustering_df %>% group_by(cluster,bac_order) %>% summarize(n=n()) %>% 
    group_by(cluster) %>% summarize(n=n()) %>% mutate(n_order = ifelse(n<4,n,"more")) %>% select(-n)
  n_phylum<-clustering_df %>% group_by(cluster,bac_phylum) %>% summarize(n=n()) %>% 
    group_by(cluster) %>% summarize(n=n()) %>% mutate(n_phylum = ifelse(n<4,n,"more")) %>% select(-n)
  n_family<-clustering_df %>% group_by(cluster,bac_family) %>% summarize(n=n()) %>% 
    group_by(cluster) %>% summarize(n=n())  %>% mutate(n_family = ifelse(n<4,n,"more")) %>% select(-n)
  n_genus<-clustering_df %>% group_by(cluster,bac_genus2) %>% summarize(n=n()) %>% 
    group_by(cluster) %>% summarize(n=n())  %>% mutate(n_genus = ifelse(n<4,n,"more")) %>% select(-n)

  n_taxa<-n_phylum %>% left_join(n_order) %>% left_join(n_family) %>% left_join(n_genus) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(-cluster,values_to = "n", names_to = "rank")
  n_taxa$rank<-factor(n_taxa$rank,levels=c("n_phylum","n_order","n_family","n_genus"))
  ggplot(n_taxa,aes(x=rank,fill=n)) %>% + geom_bar(stat="count", position ="fill")
}

### visualize: how namy clusters per taxonomic rank?
ncluster_graph<-function(clustering_df,mag_taxonomy){
  clustering_df$cluster<-clustering_df[,2]
  clustering_df <- clustering_df %>% left_join(mag_taxonomy)%>% filter(cluster>0)
  
  n_cluster_per_order<-clustering_df %>% group_by(cluster,bac_order) %>% 
    filter(bac_order!="Unknown") %>% summarize(n=n()) %>% 
    group_by(bac_order) %>% summarize(n=n()) %>% mutate(n = ifelse(n<4,n,"more")) %>% 
    mutate(rank="n_order")
  colnames(n_cluster_per_order)[1]<-"taxon"
  n_cluster_per_order$n<-as.character(n_cluster_per_order$n)
  n_cluster_per_phylum<-clustering_df %>% group_by(cluster,bac_phylum) %>% summarize(n=n()) %>% 
    group_by(bac_phylum) %>% summarize(n=n()) %>% mutate(n = ifelse(n<4,n,"more")) %>% 
    mutate(rank="n_phylum")
  colnames(n_cluster_per_phylum)[1]<-"taxon"
  n_cluster_per_phylum$n<-as.character(n_cluster_per_phylum$n)
  n_cluster_per_family<-clustering_df %>% group_by(cluster,bac_family) %>%
  filter(bac_family!="Unknown") %>% summarize(n=n()) %>% 
    group_by(bac_family) %>% summarize(n=n()) %>% mutate(n = ifelse(n<4,n,"more")) %>% 
    mutate(rank="n_family")
  colnames(n_cluster_per_family)[1]<-"taxon"
  n_cluster_per_family$n<-as.character(n_cluster_per_family$n)
  n_cluster_per_genus<-clustering_df %>% group_by(cluster,bac_genus2) %>% 
    filter(bac_genus2!="Unknown gen. sp.") %>% summarize(n=n()) %>% 
    group_by(bac_genus2) %>% summarize(n=n()) %>% mutate(n = ifelse(n<4,n,"more")) %>% 
    mutate(rank="n_genus")
  colnames(n_cluster_per_genus)[1]<-"taxon"
  n_cluster_per_genus$n<-as.character(n_cluster_per_genus$n)
  
  n_taxa<-rbind(n_cluster_per_phylum,n_cluster_per_order,n_cluster_per_family,n_cluster_per_genus) 
  n_taxa$rank<-factor(n_taxa$rank,levels=c("n_phylum","n_order","n_family","n_genus"))
  ggplot(n_taxa,aes(x=rank,fill=n)) %>% + geom_bar(stat="count", position ="fill")
}


hdbscan_graph1<-ntaxa_graph(hdbscan_results,mag_taxonomy)+ggtitle("hdbscan, # of taxa per cluster")
apcluster_graph1<-ntaxa_graph(apcluster_results,mag_taxonomy)+ggtitle("apcluster, # of taxa per cluster")
kmeans_graph1<-ntaxa_graph(kmeans_results,mag_taxonomy)+ggtitle("kmeans, # of taxa per cluster")

hdbscan_graph2<-ncluster_graph(hdbscan_results,mag_taxonomy)+ggtitle("hdbscan, # of clusters per taxon")
apcluster_graph2<-ncluster_graph(apcluster_results,mag_taxonomy)+ggtitle("apcluster, # of clusters per taxon")
kmeans_graph2<-ncluster_graph(kmeans_results,mag_taxonomy)+ggtitle("kmeans, # of clusters per taxon")

pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_taxonomic_coherence.pdf",width=8,height=8)
(hdbscan_graph1 + hdbscan_graph2) / (apcluster_graph1 + apcluster_graph2) / (kmeans_graph1 + kmeans_graph2)
dev.off()


## 10. Hierarchical clustering - based on kegg families. used all mags, including incomplete
gnm_cor = cor(t(matrix %>% select(-Genome))) 
### 10a. selecting best clustering method

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_families_clustering_comparison.pdf",width=10,height=7)
compare_clustering_methods(gnm_cor)
dev.off()

#### results: kmeans, hdbscan and apcluster look optimal

### 10b. clustering
set.seed(123)
hdbscan_gnm_fam <- cluster_by_hdbscan(gnm_cor)
hdbscan_results_fam<-data.frame(cbind(matrix$Genome,hdbscan_gnm_fam))
colnames(hdbscan_results_fam)<-c("Genome","hdbscan")

set.seed(123)
apcluster_gnm_fam <- cluster_by_apcluster(gnm_cor)
apcluster_results_fam<-data.frame(cbind(matrix$Genome,apcluster_gnm_fam))
colnames(apcluster_results_fam)<-c("Genome","apcluster")

set.seed(123)
kmeans_gnm_fam <- cluster_by_kmeans(gnm_cor)
kmeans_results_fam<-data.frame(cbind(matrix$Genome,kmeans_gnm_fam))
colnames(kmeans_results_fam)<-c("Genome","kmeans")


fam_df2<-matrix %>% left_join(mag_taxonomy) %>% 
  select(Genome,bac_family,bac_order,bac_phylum) %>%
  left_join(hdbscan_results_fam) %>%
  left_join(kmeans_results_fam) %>%
  left_join(apcluster_results_fam)
write.table(fam_df2, "analysis/07_annotate_MAGs/summarized_outputs/kegg_families_clustering.txt", sep='\t',quote = F, row.names = F, col.names = T)

##heatmap
ra_fam = rowAnnotation(df=data.frame(bac_family = fam_df2$bac_family,
                                 bac_order = fam_df2$bac_order,
                                 bac_phylum = fam_df2$bac_phylum,
                                 apcluster=fam_df2$apcluster,
                                 kmeans=fam_df2$kmeans,
                                 hdbscan=fam_df2$hdbscan)
                   #col = list(type = type_color),
                   #annotation_legend_param = list(type=list(labels=type_ordered,at=type_ordered,title = "function"))
)


row_clust = ReadDendrogram("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")
  
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_families_clustering_results.pdf",width=12,height=15)
Heatmap(data.frame(matrix) %>% select(-Genome), show_column_names = F, show_row_names = F, name = " ",
        cluster_rows = row_clust,
        right_annotation = ra_fam)
dev.off()




### 10c  taxonomic coherence
hdbscan_fam_graph1<-ntaxa_graph(hdbscan_results_fam,mag_taxonomy)+ggtitle("hdbscan, # of taxa per cluster")
apcluster_fam_graph1<-ntaxa_graph(apcluster_results_fam,mag_taxonomy)+ggtitle("apcluster, # of taxa per cluster")
kmeans_fam_graph1<-ntaxa_graph(kmeans_results_fam,mag_taxonomy)+ggtitle("kmeans, # of taxa per cluster")

hdbscan_fam_graph2<-ncluster_graph(hdbscan_results_fam,mag_taxonomy)+ggtitle("hdbscan, # of clusters per taxon")
apcluster_fam_graph2<-ncluster_graph(apcluster_results_fam,mag_taxonomy)+ggtitle("apcluster, # of clusters per taxon")
kmeans_fam_graph2<-ncluster_graph(kmeans_results_fam,mag_taxonomy)+ggtitle("kmeans, # of clusters per taxon")

pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_families_clustering_taxonomic_coherence.pdf",width=8,height=8)
(hdbscan_fam_graph1 + hdbscan_fam_graph2) / (apcluster_fam_graph1 + apcluster_fam_graph2) / (kmeans_fam_graph1 + kmeans_fam_graph2)
dev.off()

#calculate TC
l<-lapply(apcluster_results_fam$apcluster %>% unique,process_cluster,clustering_df=apcluster_results_fam,bactree=bactree)
aplcuster_fam_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(hdbscan_results_fam$hdbscan %>% unique,process_cluster,clustering_df=hdbscan_results_fam,bactree=bactree)
hdbscan_fam_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(kmeans_results_fam$kmeans %>% unique,process_cluster,clustering_df=kmeans_results_fam,bactree=bactree)
kmeans_fam_tc<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

## result: in all, the majority of clusters is not monophyletic. 
###to account for possibility of singletones, allowed exclusion of singletones
l<-lapply(apcluster_results_fam$apcluster %>% unique,process_cluster2,clustering_df=apcluster_results_fam,bactree=bactree)
aplcuster_fam_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(hdbscan_results_fam$hdbscan %>% unique,process_cluster2,clustering_df=hdbscan_results_fam,bactree=bactree)
hdbscan_fam_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)

l<-lapply(kmeans_results_fam$kmeans %>% unique,process_cluster2,clustering_df=kmeans_results_fam,bactree=bactree)
kmeans_fam_tc2<-do.call(rbind,l) %>% mutate(TC=Ngfc/Ntaxon)
##didn't improve much






e<-extract.clade(bactree,1148)



t<-apcluster_results %>% filter(apcluster==6)
mrca<-findMRCA(bactree, t$Genome, type="node")
e<-extract.clade(bactree,mrca)
plot(e)



Ancestors(bactree,t$Genome)


getMRCA(bactree,t$Genome)

listTips(bactree)[[1148]]

allDescendants(bactree,1148)







  
library(treeio)











### 9a. based on kegg families

##### 9a.1 Create Genome Functional Clusters (GFCs)
gnm_cor = cor(t(matrix %>% select(-Genome)))
gnm_cor[gnm_cor < 0] = 0 #! there are no neg correlation (as expected)



### run affinity propagation
# #!!! code to check fitting of parameter "q"
#x = seq(from=0, to=0.95, length=100)
y = sapply(x, function(i) {
   length(apcluster(s = gnm_cor, details=F, q=i,
                    lam=0.5, seed=1234, maxits=1000, convits=200)@clusters)
 })
plot(x,y, xlab="Q-vaules", ylab="# GFCs")


apcl_gnm = apcluster(s = gnm_cor, details=T, q=0.02, lam=0.5, seed=1234, maxits=1000, convits=500)
heatmap(apcl_gnm,gnm_cor)

gnm_GFC_tmp = do.call(rbind, lapply(1:length(apcl_gnm@clusters), function(i) data.frame(i, apcl_gnm@clusters[[i]])))
gnm_GFC = gnm_GFC_tmp$i[order(gnm_GFC_tmp$apcl_gnm.clusters..i.., decreasing = F)]
table(gnm_GFC)

gnm_dist = as.matrix(cophenetic(as.dendrogram(aggExCluster(s = gnm_cor, x = apcl_gnm))))
gnm_dist = gnm_dist[match(row.names(gnm_cor), row.names(gnm_dist)), 
                    match(colnames(gnm_cor), colnames(gnm_dist))]
gnm_hc = hclust(as.dist(gnm_dist), "complete")
gnm_hc = reorder(gnm_hc, wts = colSums(gnm_cor), agglo.FUN = "mean") # improve dendro sorting
gnm_hc$order = as.integer(gnm_hc$order) # otherwise it rises an issue when plotting with iheatmapr



### parse results
GFC_table = data.frame(gnm=row.names(matrix %>% select(-Genome)), GFC=gnm_GFC)
GFC_table$GFC_size = as.integer(table(GFC_table$GFC))[match(GFC_table$GFC, names(table(GFC_table$GFC)))]

###add taxonomy
GFC_table<-GFC_table %>% left_join(mag_taxonomy,by=c("gnm"="Genome"))


aggresa <- aggExCluster(gnm_cor,x=apcl_gnm)
heatmap(aggresa,gnm_cor)
plot(aggresa)







#### !!try with complete mags only
gnm_cor2 = cor(t(matrix3 %>% select(-Genome)))
x = seq(from=0, to=0.95, length=100)
y = sapply(x, function(i) {
  length(apcluster(s = gnm_cor2, details=F, q=i,
                   lam=0.5, seed=1234, maxits=1000, convits=200)@clusters)
})
plot(x,y, xlab="Q-vaules", ylab="# GFCs")


apcl_gnm2 = apcluster(s = gnm_cor2, details=T, q=0.1, lam=0.5, seed=1234, maxits=1000, convits=500)
heatmap(apcl_gnm2,gnm_cor2)



gnm_GFC_tmp = do.call(rbind, lapply(1:length(apcl_gnm2@clusters), function(i) data.frame(i, apcl_gnm2@clusters[[i]])))
gnm_GFC = gnm_GFC_tmp$i[order(gnm_GFC_tmp$apcl_gnm2.clusters..i.., decreasing = F)]
table(gnm_GFC)

gnm_dist = as.matrix(cophenetic(as.dendrogram(aggExCluster(s = gnm_cor, x = apcl_gnm))))
gnm_dist = gnm_dist[match(row.names(gnm_cor), row.names(gnm_dist)), 
                    match(colnames(gnm_cor), colnames(gnm_dist))]
gnm_hc = hclust(as.dist(gnm_dist), "complete")
gnm_hc = reorder(gnm_hc, wts = colSums(gnm_cor), agglo.FUN = "mean") # improve dendro sorting
gnm_hc$order = as.integer(gnm_hc$order) # otherwise it rises an issue when plotting with iheatmapr

### parse results
GFC_table = data.frame(gnm=row.names(matrix3 %>% select(-Genome)), GFC=gnm_GFC)
GFC_table$GFC_size = as.integer(table(GFC_table$GFC))[match(GFC_table$GFC, names(table(GFC_table$GFC)))]
###add taxonomy
GFC_table<-GFC_table %>% left_join(mag_taxonomy,by=c("gnm"="Genome"))


##### 9a.2 Create Linked Trait CLusters (LTCs)

f_r = function(x, y) {
  cont_tab = table(factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)))
  
  ## add 1 fake absence to traits without O values (M00005, M00052)
  if (sum(cont_tab[2, ]) == sum(cont_tab) | (sum(cont_tab[, 2]) == sum(cont_tab))) { # no zeros for x
    cont_tab[1, 1] = 1
  }
  ##---
  
  joint_p = cont_tab/sum(cont_tab)
  Pab = joint_p[4]
  Pa = sum(joint_p[2, ])
  Pb = sum(joint_p[, 2])
  
  D = Pab - (Pa * Pb)
  r = D / sqrt(Pa * (1 - Pa) * Pb * (1 - Pb))
  return(r)
}

pairFUN = function(mydata, fun.xy, ncore) {
  # if function ('fun.xy') generates a vector with multiple output values, specify the index ('val_index') of the desired one 
  
  smpl_comb = as.matrix(combn(x = 1:ncol(mydata), m = 2))
  
  library(future.apply)
  plan(multiprocess, workers = ncore)
  comb_out = future_apply(X = smpl_comb, MARGIN = 2, FUN = function(z) 
    fun.xy(mydata[, z[1]], mydata[, z[2]]))
  
  comb_mat = matrix(NA, nrow = ncol(mydata), ncol = ncol(mydata))
  comb_mat[lower.tri(comb_mat, diag = F)] = comb_out
  comb_mat = as.matrix(as.dist(comb_mat))
  diag(comb_mat) = 1
  dimnames(comb_mat) = list(colnames(mydata), colnames(mydata))
  
  return(comb_mat)
}

trait_cor = pairFUN(matrix3 %>% select(-Genome), fun.xy = f_r, ncore = 7)
trait_cor[trait_cor < 0] = 0

## test for significant correlation (r)
pair_r = trait_cor
pair_r[upper.tri(trait_cor, diag = T)] = NA
pair_r = melt(pair_r, na.rm = T)
pair_r$x2 = sapply(pair_r$value, function(i) i^2 * nrow(gnm_profi_filt)) # x2 test
pair_r$p.val = sapply(pair_r$x2, function(i) pchisq(i, df=2, lower.tail=F))
pair_r$p.val.adj = p.adjust(pair_r$p.val, method = "fdr")
min_signif_r = min(pair_r$value[pair_r$p.val.adj <= 0.05])
trait_cor[trait_cor < min_signif_r] = 0


##visualize
agg<-aggExCluster(gnm_cor2)
heatmap(aggres2a,gnm_cor2)
aggres2a <- aggExCluster(gnm_cor2,x=apcl_gnm2)

agg<-aggExCluster(gnm_cor2)
heatmap(aggresa,gnm_cor)
aggresa <- aggExCluster(gnm_cor,x=apcl_gnm)




Heatmap(matrix %>% select(-Genome), show_column_names = F, name = " ",)



#BiocManager::install("simplifyEnrichment",force = TRUE)








