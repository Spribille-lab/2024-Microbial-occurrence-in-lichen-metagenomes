#using multivar stats to find bacterial groups

setwd("~/Documents/coverage")


## 1. misc
library(tidyverse)
library(stringr)
library(reshape2)
library(apcluster)
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(simplifyEnrichment)
source("code/utils.R")
source("code/kegg_module_reconstruct.R")

## 2. read data
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")

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




## 5. Reconstruct kegg modules and use them for nmds. only used genomes that are >90% complete
km_dia<-KMdiagram_fetcher(create_RData=T, path="analysis/07_annotate_MAGs/rdata", new_date=Sys.Date())
km_matrix<-matrix3 %>% select(-Genome)
kmreco<-KMreco(indata=km_matrix, km_dia, len_breaks=NULL, allowed_gaps=c(1)) 
module_matrix<-kmreco[[4]]


### remove mags with low occurrence
module_df<-data.frame(module_matrix)
module_df$Genome<- matrix3$Genome
write.table( module_df, "analysis/07_annotate_MAGs/summarized_outputs/kegg_module_matrix.txt", sep='\t',quote = F, row.names = F, col.names = T)

rownames(module_df)<-module_df$Genome


## 6. Hierarchical clustering - based on kegg modules
gnm_cor_module = cor(t(module_df %>% select(-Genome))) #correlation matrix

### 6a. selecting best clustering method

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison.pdf",width=10,height=7)
compare_clustering_methods(gnm_cor_module, plot_type = "heatmap")
dev.off()

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(gnm_cor_module,plot_type = "heatmap")
dev.off()
#### results: kmeans, hdbscan and apcluster look optimal

### 6b. clustering
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


### 7c  taxonomic coherence
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


