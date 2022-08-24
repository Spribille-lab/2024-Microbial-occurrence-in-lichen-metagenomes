#cluster lichens and bacteria based on occurrence patterns

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
library(RColorBrewer)
library(phangorn)
require("phytools")
source("code/utils.R")
source("code/kegg_module_reconstruct.R")

## 2. read data
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")

## 3. get taxonomy
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_family2=ifelse(bac_family=="Unknown",paste(bac_order," fam.",sep=""),bac_family))
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family2," gen. sp.",sep=""),bac_genus))

mag_taxonomy <- mags_role %>% select(Genome,bac_genus2,bac_family2,bac_order,bac_phylum) %>% distinct()

###list of bacterial mags
bac_mags<-mags_role[!(is.na(mags_role$bac_genus)),1] %>% unique()

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont"])


## 4. make occurrence matrix
coc <- mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(Genome,metagenome,breadth) %>%
   spread(Genome, breadth, fill = 0)

M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
N = M
M[N<50] = 0#"Absent"
M[N>=50] =1# "Present"


## 5. heatmap - pre-clustering
###1. bacterial phyla, top annotation
phylum_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy,by=c("colnames.M."="Genome"))
phylum<-phylum_df$bac_phylum
ta = HeatmapAnnotation(df = data.frame(phylum = phylum),show_annotation_name=F)

###2.bacterial order, bottom annotation
order_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy,by=c("colnames.M."="Genome"))
order<-order_df$bac_order
ba = HeatmapAnnotation(df = data.frame(order = order), show_annotation_name=F)

###3. Mycobiont order. right annotation
#####prepare annotation for the mycobiont 
mtg_info$label<-mtg_info$order
lichen_group<-mtg_info %>% select(Run,label)


mycobiont_group_df<-data.frame(rownames(M)) %>% left_join(lichen_group,by=c("rownames.M."="Run"))
mycobiont_group<-mycobiont_group_df$label
myco_node_color<-c("Lecanorales" = brewer.pal(12,"Set3")[1],
                   "Peltigerales" = brewer.pal(12,"Set3")[2],
                   "Gyalectales"  = brewer.pal(12,"Set3")[3],
                   "Caliciales" = brewer.pal(12,"Set3")[4],
                   "Strigulales"  = brewer.pal(12,"Set3")[5],
                   "Verrucariales" = brewer.pal(12,"Set3")[6],
                   "Lecideales"  = brewer.pal(12,"Set3")[7],
                   "Baeomycetales"  = brewer.pal(12,"Set3")[8],
                   "Leprocaulales" = brewer.pal(12,"Set3")[9],
                   "Lichinales" = brewer.pal(12,"Set3")[10],
                   "Pertusariales" = brewer.pal(12,"Set3")[11],
                   "Ostropales" = brewer.pal(12,"Set3")[12],
                   "Umbilicariales" = brewer.pal(8,"Dark2")[1],
                   "Rhizocarpales" = brewer.pal(8,"Dark2")[2],
                   "Trypetheliales" = brewer.pal(8,"Dark2")[3],
                   "Pyrenulales"  = brewer.pal(8,"Dark2")[4],
                   "Acarosporales" = brewer.pal(8,"Dark2")[5],
                   "Teloschistales"  = brewer.pal(8,"Dark2")[6],
                   "Mycocaliciales" = brewer.pal(8,"Dark2")[7],
                   "Ostropomycetidae ins. ced." = brewer.pal(8,"Dark2")[8],
                   "Sarrameanales" = "black",
                   "Schaereriales" = "red",
                   "Arthoniales" = "green",
                   "Unknown" = brewer.pal(9,"Pastel1")[9])

###4. Sequencing depth. right annotation
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F)
depth_df<-data.frame(rownames(M)) %>% left_join(depth_df,by=c("rownames.M."="V1")) 
depth<-depth_df$V2

#### photobiont type
photobiont_df<-data.frame(rownames(M)) %>% left_join(mtg_info %>% select(Run,photobiont),by=c("rownames.M."="Run"))
photobiont<-photobiont_df$photobiont
photobiont_node_color<-c("trebouxioid" = "#4AA004",
                         "cyano" = "#4CC5C2",
                         "trebouxioid_cyano" = "#CBF9DC",
                         "absent" = "#041AA0",
                         "trentepohlioid" = "orange",
                         "Unknown" = "#CACACC",
                         "?"="#CACACC")


ra = rowAnnotation(df=data.frame(mycobiont_group = mycobiont_group, photobiont =photobiont,
                                 depth = depth),
                   col = list(mycobiont_group = myco_node_color,photobiont = photobiont_node_color),
                   annotation_legend_param = list(depth=list(title = "sequencing depth", 
                                                             at = c(0, 10000000000,20000000000,30000000000,40000000000), 
                                                             labels = c("0 bp","10 Gbp", "20 Gbp", "30 Gbp","40 Gbp")),
                                                  mycobiont_group=list(title = "mycobiont order"),photobiont=list(title = "photobiont")))

###plot
options(repr.plot.width=15, repr.plot.height=5)

HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
             bottom_annotation = ba, 
             right_annotation = ra,
             #column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))
HM


## 8. clustering: metagenomes. based on genus occurrences

matrix<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

M = as.matrix(matrix[,2:ncol(matrix)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]

###1. bacterial phyla, top annotation
phylum_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(-Genome) %>% distinct(),by=c("colnames.M."="bac_genus2"))
phylum<-phylum_df$bac_phylum
ta = HeatmapAnnotation(df = data.frame(phylum = phylum),show_annotation_name=F)

###2.bacterial order, bottom annotation
order_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(c(bac_genus2,bac_order)) %>% distinct(),by=c("colnames.M."="bac_genus2"))
order<-order_df$bac_order
ba = HeatmapAnnotation(df = data.frame(order = order), show_annotation_name=F)


HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
             bottom_annotation = ba, 
             right_annotation = ra,
             #column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))
HM


## 9. clustering: metagenomes. based on family occurrences

matrix<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_family2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_family2,values_from=occ,values_fill=0)

M = as.matrix(matrix[,2:ncol(matrix)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]

###1. bacterial phyla, top annotation
phylum_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(bac_family2,bac_phylum) %>% distinct(),by=c("colnames.M."="bac_family2"))
phylum<-phylum_df$bac_phylum
ta = HeatmapAnnotation(df = data.frame(phylum = phylum),show_annotation_name=F)

###2.bacterial order, bottom annotation
order_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(c(bac_family2,bac_order)) %>% distinct(),by=c("colnames.M."="bac_family2"))
order<-order_df$bac_order
ba = HeatmapAnnotation(df = data.frame(order = order), show_annotation_name=F)


HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
             bottom_annotation = ba, 
             right_annotation = ra,
             #column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))
HM

#### doesn't make much sense. e.g. peltigerales and cyanolichens are spread randomly. maybe the problem is with uneven depth?


## 8*. clustering: metagenomes. based on genus occurrences.
###remove rare bacteria (occurring in <5 metagenomes) and shallow metagenomes (<2Gbp)
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F)
colnames(depth_df)<-c("Run","depth")
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>% filter(depth<2000000000)

singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<5)

matrix<-mags_role_filtered %>% filter(Genome %in% bac_mags & !(bac_genus2 %in% singleton_genus$bac_genus2) & !(metagenome %in% mtg_to_exclude$Run)) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

M = as.matrix(matrix[,2:ncol(matrix)])
row_names<-matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% mutate(row_name=paste(metagenome, Lichen.metagenomes))
rownames(M) = row_names$row_name
M  =M[,colSums(M) >0 ]

###1. bacterial phyla, top annotation
phylum_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(-Genome) %>% distinct(),by=c("colnames.M."="bac_genus2"))
phylum<-phylum_df$bac_phylum
ta = HeatmapAnnotation(df = data.frame(phylum = phylum),show_annotation_name=F)

###2.bacterial order, bottom annotation
order_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(c(bac_genus2,bac_order)) %>% distinct(),by=c("colnames.M."="bac_genus2"))
order<-order_df$bac_order
ba = HeatmapAnnotation(df = data.frame(order = order), show_annotation_name=F)


#clustering
mtg_cor = cor(t(matrix %>% select(-metagenome))) #correlation matrix

set.seed(123)
pdf(file="analysis/05_MAGs/exploratory_fig/metagenomes_by_mag_occurrences_clustering_comparison.pdf",width=12,height=10)
op <- par(cex = 1.5)
compare_clustering_methods(mtg_cor,method = c("hdbscan","apcluster","MCL","walktrap","kmeans","binary_cut","dynamicTreeCut"))
dev.off()

set.seed(123)
pdf(file="analysis/05_MAGs/exploratory_fig/metagenomes_by_mag_occurrences_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(mtg_cor,method = c("hdbscan","apcluster","MCL","walktrap","kmeans","binary_cut","dynamicTreeCut"),plot_type = "heatmap")

dev.off()

bac_cor = cor(matrix %>% select(-metagenome)) #correlation matrix

set.seed(123)
pdf(file="analysis/05_MAGs/exploratory_fig/bacteria_by_mag_occurrences_clustering_comparison.pdf",width=12,height=10)
compare_clustering_methods(bac_cor,method = c("hdbscan","apcluster","MCL","walktrap","kmeans","binary_cut","dynamicTreeCut"))
dev.off()

set.seed(123)
pdf(file="analysis/05_MAGs/exploratory_fig/bacteria_by_mag_occurrences_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(bac_cor,method = c("hdbscan","apcluster","MCL","walktrap","kmeans","binary_cut","dynamicTreeCut"),plot_type = "heatmap")
dev.off()


set.seed(123)
apcluster_mtg <- cluster_by_apcluster(mtg_cor)
apcluster_mtg_results<-data.frame(cbind(matrix$metagenome,apcluster_mtg))
colnames(apcluster_mtg_results)<-c("metagenome","apcluster")

set.seed(123)
kmeans_mtg <- cluster_by_kmeans(mtg_cor)
kmeans_mtg_results<-data.frame(cbind(matrix$metagenome,kmeans_mtg))
colnames(kmeans_mtg_results)<-c("metagenome","kmeans")

mtg_clustering_df<-kmeans_mtg_results %>% left_join(apcluster_mtg_results) %>%
  left_join(mtg_info,by=c("metagenome"="Run"))
  
#write.table(module_df2, "analysis/07_annotate_MAGs/summarized_outputs/kegg_modules_clustering.txt", sep='\t',quote = F, row.names = F, col.names = T)


set.seed(123)
apcluster_bac <- cluster_by_apcluster(bac_cor)
apcluster_bac_results<-data.frame(cbind(colnames(matrix)[-1],apcluster_bac))
colnames(apcluster_bac_results)<-c("bac_genus2","apcluster")

set.seed(123)
kmeans_bac <- cluster_by_kmeans(bac_cor)
kmeans_bac_results<-data.frame(cbind(colnames(matrix)[-1],kmeans_bac))
colnames(kmeans_bac_results)<-c("bac_genus2","kmeans")

bac_clustering_df<-kmeans_bac_results %>% left_join(apcluster_bac_results) %>%
  left_join(mag_taxonomy %>% select(-Genome) %>% distinct)


##right annotation
mycobiont_group_df<-row_names %>% left_join(lichen_group,by=c("metagenome"="Run"))
mycobiont_group<-mycobiont_group_df$label
myco_node_color<-c("Lecanorales" = brewer.pal(12,"Set3")[1],
                   "Peltigerales" = brewer.pal(12,"Set3")[2],
                   "Gyalectales"  = brewer.pal(12,"Set3")[3],
                   "Caliciales" = brewer.pal(12,"Set3")[4],
                   "Strigulales"  = brewer.pal(12,"Set3")[5],
                   "Verrucariales" = brewer.pal(12,"Set3")[6],
                   "Lecideales"  = brewer.pal(12,"Set3")[7],
                   "Baeomycetales"  = brewer.pal(12,"Set3")[8],
                   "Leprocaulales" = brewer.pal(12,"Set3")[9],
                   "Lichinales" = brewer.pal(12,"Set3")[10],
                   "Pertusariales" = brewer.pal(12,"Set3")[11],
                   "Ostropales" = brewer.pal(12,"Set3")[12],
                   "Umbilicariales" = brewer.pal(8,"Dark2")[1],
                   "Rhizocarpales" = brewer.pal(8,"Dark2")[2],
                   "Trypetheliales" = brewer.pal(8,"Dark2")[3],
                   "Pyrenulales"  = brewer.pal(8,"Dark2")[4],
                   "Acarosporales" = brewer.pal(8,"Dark2")[5],
                   "Teloschistales"  = brewer.pal(8,"Dark2")[6],
                   "Mycocaliciales" = brewer.pal(8,"Dark2")[7],
                   "Ostropomycetidae ins. ced." = brewer.pal(8,"Dark2")[8],
                   "Sarrameanales" = "black",
                   "Schaereriales" = "red",
                   "Arthoniales" = "green",
                   "Unknown" = brewer.pal(9,"Pastel1")[9])

###4. Sequencing depth. right annotation
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F)
depth_df<-mtg_clustering_df %>% left_join(depth_df,by=c("metagenome"="V1")) 
depth<-depth_df$V2

#### photobiont type
photobiont<-row_names$photobiont
photobiont_node_color<-c("trebouxioid" = "#4AA004",
                         "cyano" = "#4CC5C2",
                         "trebouxioid_cyano" = "#CBF9DC",
                         "absent" = "#041AA0",
                         "trentepohlioid" = "orange",
                         "Unknown" = "#CACACC",
                         "?"="#CACACC")


ra2 = rowAnnotation(df=data.frame(mycobiont_group = mtg_clustering_df$order, photobiont =mtg_clustering_df$photobiont,
      depth = depth,apcluster=mtg_clustering_df$apcluster,kmeans=mtg_clustering_df$kmeans),
                   col = list(mycobiont_group = myco_node_color,photobiont = photobiont_node_color),
                   annotation_legend_param = list(depth=list(title = "sequencing depth", 
                                                             at = c(0, 10000000000,20000000000,30000000000,40000000000), 
                                                             labels = c("0 bp","10 Gbp", "20 Gbp", "30 Gbp","40 Gbp"))))

###1. bacterial phyla, top annotation
phylum_df<-data.frame(colnames(M)) %>% left_join(mag_taxonomy %>% select(bac_genus2,bac_phylum) %>% distinct(),by=c("colnames.M."="bac_genus2"))
phylum<-phylum_df$bac_phylum
ta = HeatmapAnnotation(df = data.frame(phylum = phylum),show_annotation_name=F)

###2.bacterial clusters, bottom annotation
bac_clustering_df2 <-data.frame(colnames(M)) %>% left_join(bac_clustering_df %>% select(c(bac_genus2,kmeans, apcluster)) %>% distinct(),by=c("colnames.M."="bac_genus2"))

ba = HeatmapAnnotation(df = data.frame(kmeans = bac_clustering_df2$kmeans,apcluster=bac_clustering_df2$apcluster), show_annotation_name=F)


HM = Heatmap(M, show_row_names = T, show_column_names = T, name = " ",
             top_annotation = ta, 
             bottom_annotation = ba, 
             right_annotation = ra2,
             #clustering_method_columns = "apcluster", 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))

pdf(file="analysis/05_MAGs/exploratory_fig/clustering_metagenomes_vs_bac_genera.pdf",width=15,height=30)
HM
dev.off()
 ##this is the version of clustering and heatmap that more or less worked!




## 10 CCA, based on families of bacteria. Used phtoboint type, mycobiont order, and depth as predictors
matrix<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_family2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_family2,values_from=occ,values_fill=0)

### 10a all metagenomes
library(vegan)
mtg_to_exclude<-mtg_info %>% filter(photobiont=="Unknown" | order =="Unknown" | is.na(order))

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))

matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

summary(cca_fam)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

### 10b remove metagenomes with <2Gbp
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F,col.names = c("Run","depth"))
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order =="Unknown" | is.na(order) | depth<2000000000)

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

summary(cca_fam)
plot(cca_fam)
plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")


## 10d Lichinales and Verrucariales are outliers. tried to remove them
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order %in% c("Unknown","Lichinales","Verrucariales") | is.na(order) | depth<2000000000)

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

##10e now Umbilicariales is an outlier. removed it too
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order %in% c("Unknown","Lichinales","Verrucariales","Umbilicariales") | is.na(order) | depth<2000000000)

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

## 10f remove bacterial families that occurred in only one metagenome. remove all metagenomes <2 Gbp
singleton_fam<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(metagenome,bac_family2) %>% distinct() %>%
  group_by(bac_family2) %>% summarize(n=n()) %>% filter(n==1)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order =="Unknown" | is.na(order) | depth<2000000000)

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

## 10g Verrucariales and Licheinales still are outliers. removed them
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order  %in% c("Unknown","Verrucariales","Lichinales") | is.na(order) | depth<2000000000)

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

pdf(file="analysis/05_MAGs/exploratory_fig/cca_metagenomes_vs_bac_families_removed_outliers.pdf",width=7,height=7)
plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2, removed Verrucariales and Lichinales")
dev.off()


### 10h. check just peltigerales
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order != "Peltigerales" | is.na(order) | depth<2000000000)

matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2))
not_occ_fam<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2),-one_of(not_occ_fam))

matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2),-one_of(not_occ_fam)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_fam<-cca(matrix_for_cca[,-1] ~ + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

##10i. remove bacterial families that occurred in <5 metagenomes. remove all metagenomes <2 Gbp
singleton_fam<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(metagenome,bac_family2) %>% distinct() %>%
  group_by(bac_family2) %>% summarize(n=n()) %>% filter(n<5)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order =="Unknown" | is.na(order) | depth<2000000000)

matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2)) %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% filter(sum>0) %>% select(-sum)
matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(metagenome %in% matrix_for_cca$metagenome) %>% select(-one_of(singleton_fam$bac_family2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")

#remove lichinales
mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order  %in% c("Unknown","Lichinales") | is.na(order) | depth<2000000000)

matrix_for_cca<-matrix %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_fam$bac_family2)) %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% filter(sum>0) %>% select(-sum)
matrix_df <- matrix %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(metagenome %in% matrix_for_cca$metagenome) %>% select(-one_of(singleton_fam$bac_family2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_fam<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_fam, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2")



## 11. CCA based on bacterial genera

### 11a. on all metagenoems with >2 Gbp of data
matrix_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order =="Unknown" | is.na(order) | depth<2000000000)

### 11b. excludegenera that occurred <2 times
singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<2)

matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2))
not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus))

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_genus<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

library(ggvegan)
autoplot(cca_genus)

### 11c Verrucariales is an outlier

matrix_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order  %in% c("Unknown","Verrucariales","Lichinales") | is.na(order) | depth<2000000000)

#excludegenera that occurred <2 times
singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<2)

matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2))
not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus))

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_genus<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)


autoplot(cca_genus)



### 11d. check peltigerales by bacterial genus
matrix_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order != "Peltigerales" | depth<2000000000)
peltigerales_list<-mtg_info %>% filter(order=="Peltigerales")


### 11e. excludegenera that occurred <2 times in the peltigerales subset
singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  filter(metagenome %in% peltigerales_list$Run) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<2)

matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2))
not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus))

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_genus<-cca(matrix_for_cca[,-1] ~ + photobiont + depth, data = matrix_df)

pdf(file="analysis/05_MAGs/exploratory_fig/cca_peltigerales_vs_bac_genera.pdf",width=7,height=7)
plot(cca_genus, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2. peltigerales")
dev.off()


anova(cca_genus, by="axis")
#results:
#Df ChiSquare      F Pr(>F)
#CCA1      1    0.2076 0.9999  0.716
#CCA2      1    0.1808 0.8708  0.650
#Residual 16    3.3221              


### 11e. excludegenera that occurred <4 times in the peltigerales subset
singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  filter(metagenome %in% peltigerales_list$Run) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<4)

matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2))
not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus))

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_genus<-cca(matrix_for_cca[,-1] ~ + photobiont + depth, data = matrix_df)

pdf(file="analysis/05_MAGs/exploratory_fig/cca_peltigerales_vs_bac_genera.pdf",width=7,height=7)
plot(cca_genus, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2. peltigerales")
dev.off()


anova(cca_genus, by="axis")
#results:
#Df ChiSquare      F Pr(>F)
#CCA1      1   0.13036 1.3072  0.416
#CCA2      1   0.02594 0.2602  0.998
#Residual 16   1.59555              



### 11b. exclude genera that occurred <5 times
matrix_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order =="Unknown" | is.na(order) | depth<2000000000)

singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  filter(metagenome %in% peltigerales_list$Run) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<5)

not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% 
  select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>% filter(sum>0) %>% select(-sum)

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(metagenome %in% matrix_for_cca$metagenome) %>% select(-one_of(singleton_genus$bac_genus2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))

cca_genus<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_genus, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2. remove rare genera")


##remove outliers lichinales and umbilicariales
matrix_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>% 
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  mutate(occ = ifelse(n>0,1,0)) %>%
  pivot_wider(-n,names_from=bac_genus2,values_from=occ,values_fill=0)

mtg_to_exclude<-mtg_info %>% left_join(depth_df) %>%
  filter(photobiont=="Unknown" | order %in% c("Unknown","Lichinales", "Umbilicariales","Verrucariales","Acarosporales" ) | is.na(order) | depth<2000000000)

singleton_genus<-mags_role_filtered %>% filter(Genome %in% bac_mags) %>%
  filter(metagenome %in% peltigerales_list$Run) %>%
  select(metagenome,bac_genus2) %>% distinct() %>%
  group_by(bac_genus2) %>% summarize(n=n()) %>% filter(n<5)

not_occ_genus<-colnames(matrix_for_cca %>% select(-metagenome))[colSums(matrix_for_cca %>% select(-metagenome)) == 0]
matrix_for_cca<-matrix_genus %>% filter(!(metagenome %in% mtg_to_exclude$Run)) %>% 
  select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>% filter(sum>0) %>% select(-sum)

matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(metagenome %in% matrix_for_cca$metagenome) %>% select(-one_of(singleton_genus$bac_genus2)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))

cca_genus<-cca(matrix_for_cca[,-1] ~ order + photobiont + depth, data = matrix_df)

plot(cca_genus, scaling = 2, 
     display = c("species", "cn"), 
     main = "biplot cca, scaling 2. remove rare genera, remove lichinales and umbilicariales")










matrix_df <- matrix_genus %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% 
  filter(!(metagenome %in% mtg_to_exclude$Run)) %>% select(-one_of(singleton_genus$bac_genus2),-one_of(not_occ_genus)) %>%
  left_join(depth_df,by=c("metagenome"="Run"))
cca_genus<-cca(matrix_for_cca[,-1] ~ + photobiont + depth, data = matrix_df)












## 6. clustering: metagenomes. based on MAG occurrences
mtg_cor = cor(t(M)) #correlation matrix

### 6a. selecting best clustering method

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison.pdf",width=10,height=7)
compare_clustering_methods(mtg_cor)
dev.off()

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(mtg_cor,plot_type = "heatmap")
dev.off()
#### results: kmeans, hdbscan and apcluster look optimal

## 7. clustering: metagenomes. based on MAG occurrences
mtg_cor = cor(t(M)) #correlation matrix

### 7a. selecting best clustering method

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison.pdf",width=10,height=7)
compare_clustering_methods(mtg_cor)
dev.off()

set.seed(123)
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison2.pdf",width=10,height=7)
compare_clustering_methods(mtg_cor,plot_type = "heatmap")
dev.off()
#### results: kmeans, hdbscan and apcluster look optimal

