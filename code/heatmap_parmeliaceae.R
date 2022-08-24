setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(RColorBrewer)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
depth_df<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F)

#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$bac_phylum=="Cyanobacteria"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

mags_role$label[mags_role$bac_phylum=="Actinobacteriota"]<-"Actinobacteriota"

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

# filter: select only parmeliaceae and metagenomes with >2gbp

mtg_info<-mtg_info%>% left_join(depth_df,by=c("Run"="V1"))
selected<-mtg_info %>% filter(V2>=2000000000, family=="Parmeliaceae")

mags_role_filtered<-mags_role_filtered %>% filter(metagenome %in% selected$Run)


#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")
#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- ape::drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722173_concoct_merged.0_mycobiont_Melanohalea_halei_Melanohalea_halei"),resolve.root = TRUE)



#prepped tree for the heatmap
ta = table(mycobiont_mags)
mtg_to_show<-mags_role_filtered %>% filter(confirmed_role=="mycobiont")

modded_tree = mycobiont_tree
#changed tip names to match the Genome column
n1<-str_split_fixed(modded_tree$tip.label, "_", 5) %>% data.frame() %>% mutate(name=paste(X1,X2,X3,X4,sep="_"))
modded_tree$tip.label <- n1$name

for (species in names(ta)){
  leafes_needed  = sum(mtg_to_show$Genome == species)
  while(sum(modded_tree$tip.label == species) < leafes_needed){
    modded_tree = suppressWarnings(bind.tip(modded_tree, species, edge.length=0.0001, where=which(modded_tree$tip.label==species)[1], position=0))
  }
  modded_tree$tip.label[modded_tree$tip.label == species] = mtg_to_show$metagenome[mtg_to_show$Genome == species]
  
}
write.tree(modded_tree, "analysis/05_MAGs/tmptree.tre")
row_clust = ReadDendrogram("analysis/05_MAGs/tmptree.tre")

##load bacterial tree
bactree <- read.tree("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")


#function to produce heatmap for bacterial groups
make_hm_bac<-function(mags_role_filtered,group_to_include,bactree,row_clust,node_color){
  
  #made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
  coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label) %>%
    filter(label %in% group_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
    mutate(label = NULL) %>% spread(Genome, breadth, fill = 0) %>% 
    select(-one_of(mycobiont_mags$Genome)) #remove mycobiont mags
  
  
  #drop unused tips
  tip_to_remove<-bactree$tip.label[!(bactree$tip.label %in% colnames(coc))]
  
  bactree_filtered<-ape::drop.tip(bactree, tip_to_remove)
  write.tree(bactree_filtered, "analysis/05_MAGs/tmptree2.tre")
  col_clust = ReadDendrogram("analysis/05_MAGs/tmptree2.tre")
  
  #make the matrix
  M = as.matrix(coc[,2:ncol(coc)])
  coc2<-coc %>% left_join(mtg_info,by=c("metagenome"="Run")) %>% mutate(row_label=paste(metagenome,Lichen.metagenomes,sep=" "))
  
  rownames(M) = coc2$row_label
  M  =M[,colSums(M) >0 ]
  N = M
  M[N<50] = 0#"Absent"
  M[N>=50] =1# "Present"
  
  o1 = seriate(dist(M), method = "TSP")
  o2 = seriate(dist(t(M)), method = "TSP")
  
  #annotations:
  #1. group, top annotation
  clade = rep(group_to_include, ncol(M))
  ta = HeatmapAnnotation(df = data.frame(group = clade), col = list(group = node_color),show_annotation_name=F)
  
  #2.bacterial genus, bottom annotation
  genus_df<-data.frame(colnames(M)) %>% left_join(mags_role_filtered %>%select(Genome, bac_genus) %>%unique(),by=c("colnames.M."="Genome"))
  genus<-genus_df$bac_genus
  ba = HeatmapAnnotation(df = data.frame(genus = genus), show_annotation_name=F)
  
 
  #4. Sequencing depth. right annotation
    #depth_df2<-coc2 %>% left_join(depth_df,by=c("metagenome"="V1")) 
  depth<-coc2$V2
  
  ra = rowAnnotation(df=data.frame(depth = depth),
                     annotation_legend_param = list(depth=list(title = "sequencing depth", 
                                                               at = c(0, 10000000000,20000000000,30000000000,40000000000), 
                                                               labels = c("0 bp","10 Gbp", "20 Gbp", "30 Gbp","40 Gbp"))
                                                   ))
  
  #plot
  options(repr.plot.width=15, repr.plot.height=5)
  
  HM = Heatmap(M, show_row_names = T, show_column_names = F, name = " ",
               top_annotation = ta, 
               bottom_annotation = ba, 
               right_annotation = ra,
               column_order = get_order(o2),
               #column_split = clade, 
               #cluster_rows = row_clust,
               rect_gp = gpar(col= "#cecece"), 
               cluster_rows = reorder(row_clust, get_order(o1)) ,
               cluster_columns = reorder(col_clust, get_order(o2)) ,
               row_dend_width = unit(8, "cm"),
               #cell_fun = function(j, i, x, y, width, height, fill) {
               #          grid.rect(x = x, y = y, width = width, height = height, 
               #            gp = gpar(col = 'white', fill = 'white'))
               #          grid.circle(x = x, y = y, r = 0.005,
               #                gp = gpar(fill = fill, col = "#cecece"))
               #    }    ,
               col = c("0" = "#ffffff", "1" = "#000000"))
  
  
  return(HM)
}


#produce heatmaps

#groups_to_include<-c("Cyanobacteria","Green Algal Photobiont","Acidobacteriaceae","Beijerinckiaceae","Acetobacteraceae")
group_to_include<-c("Acetobacteraceae")
node_color<-c("Acetobacteraceae" = "#c43b0e")
hm_aceto<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)

pdf(file="analysis/05_MAGs/exploratory_fig/parmeliaceae_aceto.pdf",width=10,height=7)
hm_aceto
dev.off()


group_to_include<-c("Acidobacteriaceae")
node_color<-c("Acidobacteriaceae" = "#6c44da")
hm_acido<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)

pdf(file="analysis/05_MAGs/exploratory_fig/parmeliaceae_acido.pdf",width=10,height=7)
hm_acido
dev.off()


group_to_include<-c("Beijerinckiaceae")
node_color<-c("Beijerinckiaceae" = "#783f04")
hm_beij<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)

pdf(file="analysis/05_MAGs/exploratory_fig/parmeliaceae_beij.pdf",width=10,height=7)
hm_beij
dev.off()


pdf(file="analysis/05_MAGs/exploratory_fig/parmeliaceae_beij_aceto.pdf",width=10,height=7)
hm_beij+hm_aceto
dev.off()


group_to_include<-c("UBA10450")
node_color<-c("UBA10450" = "#1B9E77")
hm_UBA10450<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)
pdf(file="analysis/05_MAGs/exploratory_fig/parmeliaceae_UBA10450.pdf",width=10,height=7)
hm_UBA10450
dev.off()



group_to_include<-c("Cyanobacteria")
node_color<-c("Cyanobacteria" = "#09b1db")
hm_cyano<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)


group_to_include<-c("Green Algal Photobiont")
node_color<-c("Green Algal Photobiont" = "#00CD6C")
euktree_intact_names<-read.tree("analysis/05_MAGs/trees/eukaryotes.treefile")
euktree_intact_names$tip.label<-str_replace(euktree_intact_names$tip.label, ".fa", "")
hm_alga<-make_hm_bac(mags_role_filtered,group_to_include,bactree=euktree_intact_names,row_clust,node_color)


group_to_include<-c("Cyphobasidium")
node_color<-c("Cyphobasidium"="#c90076")
euktree_intact_names<-read.tree("analysis/05_MAGs/trees/eukaryotes.treefile")
euktree_intact_names$tip.label<-str_replace(euktree_intact_names$tip.label, ".fa", "")
hm_cypho<-make_hm_bac(mags_role_filtered,group_to_include,bactree=euktree_intact_names,row_clust,node_color)

group_to_include<-c("Tremella")
node_color<-c("Tremella"="#e091bf")
euktree_intact_names<-read.tree("analysis/05_MAGs/trees/eukaryotes.treefile")
euktree_intact_names$tip.label<-str_replace(euktree_intact_names$tip.label, ".fa", "")
hm_trem<-make_hm_bac(mags_role_filtered,group_to_include,bactree=euktree_intact_names,row_clust,node_color)
hm_alga+hm_trem

group_to_include<-c("Actinobacteriota")
node_color<-c("Actinobacteriota" = "#c43b0e")
hm_actino<-make_hm_bac(mags_role_filtered,group_to_include,bactree,row_clust,node_color)


hm_cyano+hm_sph

hm_beij
hm_aceto
hm_acido
hm_alga
