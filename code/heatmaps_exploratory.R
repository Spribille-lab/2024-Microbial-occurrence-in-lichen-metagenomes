#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

#groups to include in the heatmap
groups_to_include<-c("Cyanobacteria","Green Algal Photobiont","Acidobacteriaceae","Beijerinckiaceae","Acetobacteraceae")

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")


#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722280_metabat2_bin.4_mycobiont_Dictyomeridium_proponens",
"public_SRR14721933_metabat2_bin.1_mycobiont_Bathelium_carolinianum",
"public_SRR14722281_metabat2_bin.1_mycobiont_Viridothelium_virens",
"public_SRR14722143_metabat2_bin.4_mycobiont_Zwackhia_viridis",
"public_SRR14721966_metabat2_bin.6_mycobiont_Opegrapha_moroziana",
"public_SRR14721982_metabat2_bin.3_mycobiont_fungal_sp._Lendemer_49042_NY-3033146",
"public_SRR14722233_metabat2_bin.15_mycobiont_Arthothelium_ruanum",
"public_SRR14721978_concoct_merged.0_mycobiont_Arthonia_rubella",
"public_SRR13685140_metabat2_bin.1_mycobiont_Arthonia_susa"),resolve.root = TRUE)



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

#made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
  mutate(label = NULL) %>% spread(Genome, breadth, fill = 0) %>% 
  select(-one_of(mycobiont_mags$Genome)) #remove mycobiont mags

#made the matrix
M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
N = M
M[N<50] = 0#"Absent"
M[N>=50] =1# "Present"
clade = rep(NA, ncol(M))
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Cyanobacteria"))$Genome] = "Cyanobacteria"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Green Algal Photobiont"))$Genome] = "Green Algal Photobiont"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Acidobacteriaceae"))$Genome] = "Acidobacteriaceae"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Beijerinckiaceae"))$Genome] = "Beijerinckiaceae"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Acetobacteraceae"))$Genome] = "Acetobacteraceae"
o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=15, repr.plot.height=5)
node_colors = c("Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria" = "#09b1db",                  
                "Fungi" = "#FFC61E",
                "Acetobacteraceae" = "#c43b0e",
                "Beijerinckiaceae" = "#783f04",
                "Acidobacteriaceae" = "#6c44da"
)

ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
             column_order = get_order(o2),
             column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             cluster_rows = reorder(row_clust, get_order(o1)) ,
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))


HM




###heatmap for selected bacterial genera
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$label<-mags_role$bac_genus
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

#groups to include in the heatmap
groups_to_include<-c("Cyanobacteria","Green Algal Photobiont","Lichenihabitans","CAHJXG01","LMUY01","EB88","RH-AL1")

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")


#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722280_metabat2_bin.4_mycobiont_Dictyomeridium_proponens",
                                                 "public_SRR14721933_metabat2_bin.1_mycobiont_Bathelium_carolinianum",
                                                 "public_SRR14722281_metabat2_bin.1_mycobiont_Viridothelium_virens",
                                                 "public_SRR14722143_metabat2_bin.4_mycobiont_Zwackhia_viridis",
                                                 "public_SRR14721966_metabat2_bin.6_mycobiont_Opegrapha_moroziana",
                                                 "public_SRR14721982_metabat2_bin.3_mycobiont_fungal_sp._Lendemer_49042_NY-3033146",
                                                 "public_SRR14722233_metabat2_bin.15_mycobiont_Arthothelium_ruanum",
                                                 "public_SRR14721978_concoct_merged.0_mycobiont_Arthonia_rubella",
                                                 "public_SRR13685140_metabat2_bin.1_mycobiont_Arthonia_susa"),resolve.root = TRUE)



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

#made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
  mutate(label = NULL) %>% spread(Genome, breadth, fill = 0) %>% 
  select(-one_of(mycobiont_mags$Genome)) #remove mycobiont mags

#made the matrix
M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
N = M
M[N<50] = 0#"Absent"
M[N>=50] =1# "Present"
clade = rep(NA, ncol(M))
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Cyanobacteria"))$Genome] = "Cyanobacteria"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Green Algal Photobiont"))$Genome] = "Green Algal Photobiont"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Lichenihabitans"))$Genome] = "Lichenihabitans"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="CAHJXG01"))$Genome] = "CAHJXG01"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="LMUY01"))$Genome] = "LMUY01"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="EB88"))$Genome] = "EB88"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="RH-AL1"))$Genome] = "RH-AL1"

o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=15, repr.plot.height=5)
node_colors = c("Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria" = "#09b1db",                  
                "Fungi" = "#FFC61E",
                "Lichenihabitans" = "#c43b0e",
                "EB88" = "#783f04",
                "CAHJXG01" = "#6c44da",
                "LMUY01" = "#FFC61E",
                "RH-AL1" = "#969603"
)

ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
             column_order = get_order(o2),
             column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             cluster_rows = reorder(row_clust, get_order(o1)) ,
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000"))

png("results/figures/heatmap_full.png")
HM

###heatmap on the bacterial genus level
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(RColorBrewer)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_fam <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")


mags_role$label<-mags_role$bac_fam
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

#groups to include in the heatmap
groups_to_include<-c("Sphingomonadaceae","Acidobacteriaceae","Beijerinckiaceae","Acetobacteraceae")

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")


#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722280_metabat2_bin.4_mycobiont_Dictyomeridium_proponens",
                                                 "public_SRR14721933_metabat2_bin.1_mycobiont_Bathelium_carolinianum",
                                                 "public_SRR14722281_metabat2_bin.1_mycobiont_Viridothelium_virens",
                                                 "public_SRR14722143_metabat2_bin.4_mycobiont_Zwackhia_viridis",
                                                 "public_SRR14721966_metabat2_bin.6_mycobiont_Opegrapha_moroziana",
                                                 "public_SRR14721982_metabat2_bin.3_mycobiont_fungal_sp._Lendemer_49042_NY-3033146",
                                                 "public_SRR14722233_metabat2_bin.15_mycobiont_Arthothelium_ruanum",
                                                 "public_SRR14721978_concoct_merged.0_mycobiont_Arthonia_rubella",
                                                 "public_SRR13685140_metabat2_bin.1_mycobiont_Arthonia_susa"),resolve.root = TRUE)



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

#made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label,bac_genus) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
 group_by(metagenome,bac_genus) %>% summarize(n=n()) %>%
 spread(bac_genus, n, fill = 0) %>% select(-c(Unknown,"<NA>"))




#made the matrix
M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
#N = M
#M[N<50] = 0#"Absent"
#M[N>=50] =1# "Present"

#get table where each genus is matched to its family
classification<-mags_role_filtered %>% select(bac_genus,bac_fam) %>% unique()
clade = rep(NA, ncol(M))
clade[colnames(M)  %in% (classification %>% filter(bac_fam=="Sphingomonadaceae"))$bac_genus] = "Sphingomonadaceae"
clade[colnames(M)  %in% (classification %>% filter(bac_fam=="Acidobacteriaceae"))$bac_genus] = "Acidobacteriaceae"
clade[colnames(M)  %in% (classification %>% filter(bac_fam=="Beijerinckiaceae"))$bac_genus] = "Beijerinckiaceae"
clade[colnames(M)  %in% (classification %>% filter(bac_fam=="Acetobacteraceae"))$bac_genus] = "Acetobacteraceae"
o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=30, repr.plot.height=10)
node_colors = c("Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria" = "#09b1db",                  
                "Fungi" = "#FFC61E",
                "Acetobacteraceae" = "#c43b0e",
                "Beijerinckiaceae" = "#783f04",
                "Acidobacteriaceae" = "#6c44da",
                "Sphingomonadaceae" = "#34b7eb"
)

ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
cols<-colorRampPalette(brewer.pal(9, "Reds"))(9)
HM = Heatmap(M, show_row_names = F, show_column_names = T, name = " ",
             top_annotation = ta, 
             column_order = get_order(o2),
             column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             cluster_rows = reorder(row_clust, get_order(o1)) ,
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("white",colorRampPalette(brewer.pal(9, "Reds"))(9)[3:9],"#3b0108")
             )
HM
png("analysis/05_MAGs/exploratory_fig/heatmap_bac_genus_level.png")

####add mag trees for the columns
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

#groups to include in the heatmap
#groups_to_include<-c("Cyanobacteria","Green Algal Photobiont","Acidobacteriaceae","Beijerinckiaceae","Acetobacteraceae")
groups_to_include<-c("Cyanobacteria","Acidobacteriaceae","Beijerinckiaceae","Acetobacteraceae")


#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")


#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722280_metabat2_bin.4_mycobiont_Dictyomeridium_proponens",
                                                 "public_SRR14721933_metabat2_bin.1_mycobiont_Bathelium_carolinianum",
                                                 "public_SRR14722281_metabat2_bin.1_mycobiont_Viridothelium_virens",
                                                 "public_SRR14722143_metabat2_bin.4_mycobiont_Zwackhia_viridis",
                                                 "public_SRR14721966_metabat2_bin.6_mycobiont_Opegrapha_moroziana",
                                                 "public_SRR14721982_metabat2_bin.3_mycobiont_fungal_sp._Lendemer_49042_NY-3033146",
                                                 "public_SRR14722233_metabat2_bin.15_mycobiont_Arthothelium_ruanum",
                                                 "public_SRR14721978_concoct_merged.0_mycobiont_Arthonia_rubella",
                                                 "public_SRR13685140_metabat2_bin.1_mycobiont_Arthonia_susa"),resolve.root = TRUE)



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


#made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
  mutate(label = NULL) %>% spread(Genome, breadth, fill = 0) %>% 
  select(-one_of(mycobiont_mags$Genome)) #remove mycobiont mags

##load bacterial tree
bactree <- read.tree("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")
#drop unused tips
tip_to_remove<-bactree$tip.label[!(bactree$tip.label %in% colnames(coc))]

bactree_filtered<-drop.tip(bactree, tip_to_remove)
write.tree(bactree_filtered, "analysis/05_MAGs/tmptree2.tre")
col_clust = ReadDendrogram("analysis/05_MAGs/tmptree2.tre")

#made the matrix
M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
N = M
M[N<50] = 0#"Absent"
M[N>=50] =1# "Present"
clade = rep(NA, ncol(M))
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Cyanobacteria"))$Genome] = "Cyanobacteria"
#clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Green Algal Photobiont"))$Genome] = "Green Algal Photobiont"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Acidobacteriaceae"))$Genome] = "Acidobacteriaceae"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Beijerinckiaceae"))$Genome] = "Beijerinckiaceae"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Acetobacteraceae"))$Genome] = "Acetobacteraceae"
o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=15, repr.plot.height=5)
node_colors = c("Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria" = "#09b1db",                  
                "Fungi" = "#FFC61E",
                "Acetobacteraceae" = "#c43b0e",
                "Beijerinckiaceae" = "#783f04",
                "Acidobacteriaceae" = "#6c44da"
)

ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
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


HM

###heatmap for selected bacterial genera with tree for the columns
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))

#load the data set
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$label<-mags_role$bac_genus
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

#groups to include in the heatmap
groups_to_include<-c("Cyanobacteria","Lichenihabitans","CAHJXG01","LMUY01","EB88","RH-AL1")

#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])

#load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile")


#drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% filter(confirmed_role=="mycobiont") %>% select(Genome) %>%
  unique()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14722280_metabat2_bin.4_mycobiont_Dictyomeridium_proponens",
                                                 "public_SRR14721933_metabat2_bin.1_mycobiont_Bathelium_carolinianum",
                                                 "public_SRR14722281_metabat2_bin.1_mycobiont_Viridothelium_virens",
                                                 "public_SRR14722143_metabat2_bin.4_mycobiont_Zwackhia_viridis",
                                                 "public_SRR14721966_metabat2_bin.6_mycobiont_Opegrapha_moroziana",
                                                 "public_SRR14721982_metabat2_bin.3_mycobiont_fungal_sp._Lendemer_49042_NY-3033146",
                                                 "public_SRR14722233_metabat2_bin.15_mycobiont_Arthothelium_ruanum",
                                                 "public_SRR14721978_concoct_merged.0_mycobiont_Arthonia_rubella",
                                                 "public_SRR13685140_metabat2_bin.1_mycobiont_Arthonia_susa"),resolve.root = TRUE)



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

#made the occurrence matrix, showing all metgaenomes which yielded the mycobiont MAG and all Genomes from the groups of interest
coc <- mags_role_filtered %>% select(Genome,metagenome,breadth,label) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% filter(metagenome %in% mtg_to_show$metagenome) %>%
  mutate(label = NULL) %>% spread(Genome, breadth, fill = 0) %>% 
  select(-one_of(mycobiont_mags$Genome)) #remove mycobiont mags


##load bacterial tree
bactree <- read.tree("analysis/05_MAGs/trees/gtdbtk.bac120.user_msa.fasta.treefile")
#drop unused tips
tip_to_remove<-bactree$tip.label[!(bactree$tip.label %in% colnames(coc))]

bactree_filtered<-drop.tip(bactree, tip_to_remove)
write.tree(bactree_filtered, "analysis/05_MAGs/tmptree2.tre")
col_clust = ReadDendrogram("analysis/05_MAGs/tmptree2.tre")


#made the matrix
M = as.matrix(coc[,2:ncol(coc)])
rownames(M) = coc$metagenome
M  =M[,colSums(M) >0 ]
N = M
M[N<50] = 0#"Absent"
M[N>=50] =1# "Present"
clade = rep(NA, ncol(M))
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Cyanobacteria"))$Genome] = "Cyanobacteria"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Green Algal Photobiont"))$Genome] = "Green Algal Photobiont"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="Lichenihabitans"))$Genome] = "Lichenihabitans"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="CAHJXG01"))$Genome] = "CAHJXG01"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="LMUY01"))$Genome] = "LMUY01"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="EB88"))$Genome] = "EB88"
clade[colnames(M)  %in% (mags_role_filtered %>% filter(label=="RH-AL1"))$Genome] = "RH-AL1"

o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=15, repr.plot.height=5)
node_colors = c("Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria" = "#09b1db",                  
                "Fungi" = "#FFC61E",
                "Lichenihabitans" = "#c43b0e",
                "EB88" = "#783f04",
                "CAHJXG01" = "#6c44da",
                "LMUY01" = "#FFC61E",
                "RH-AL1" = "#969603"
)

ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
HM = Heatmap(M, show_row_names = F, show_column_names = F, name = " ",
             top_annotation = ta, 
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

png("results/figures/heatmap_full.png")
HM
