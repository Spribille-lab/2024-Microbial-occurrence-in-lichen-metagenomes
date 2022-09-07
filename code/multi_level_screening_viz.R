# Heatmap showing detection of symbionts by mag, assembly, reads
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
library(conflicted)
options(repr.plot.width=20, repr.plot.height=20)
theme_set(theme_minimal(base_size = 23))
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# 1. Make a table for presence/absence of mags

#load the mag info
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Trebouxiophyceae"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cystobasidiomycetes"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremellomycetes"

#groups to include in the heatmap
groups_to_include<-c("Trebouxiophyceae","Acidobacteriaceae",
                     "Beijerinckiaceae","Acetobacteraceae",
                     "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes")

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

mycobiont_tree <- ape::drop.tip(tree, tree$tip.label[-ind[,2]])
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


mag_presence <- mags_role_filtered %>% select(metagenome,breadth,label) %>%
  filter(label %in% groups_to_include | label=="Mycobiont") %>% 
  filter(metagenome %in% mtg_to_show$metagenome) %>%
  filter(breadth>50) %>% select(-breadth) %>% distinct() %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label) %>%
  select(-Mycobiont) 

# 2. get info from metaxa for 18S of CYpho, Tremella, and Trebouxia. screening of assemblies and reads
matrix_level5<-read.delim("analysis/03_metagenome_reanalysis/metaxa_level_5_combined.txt")
matrix_level5_long<-matrix_level5 %>% pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% mutate(occurrence=ifelse(abundance>0,1,0))

#summarize by group, reduce the table only to target groups 
metaxa_occ_target<-matrix_level5_long %>% filter(occurrence>0, grepl( "Cystobasidiomycetes", Taxa) |
                                grepl( "Tremellomycetes", Taxa) |
                                grepl( "Trebouxi", Taxa)) %>%
  mutate(label=ifelse(grepl("Cystobasidiomycetes",Taxa),"Cystobasidiomycetes",
               ifelse(grepl("Beijerinckiaceae",Taxa), "Beijerinckiaceae",
               ifelse(grepl("Tremellomycetes", Taxa),"Tremellomycetes",
               ifelse(grepl("Trebouxi", Taxa),"Trebouxiophyceae",NA))))) %>%
  select(-Taxa,-abundance) %>% group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))


## 3. get info from idtaxa for 16S screening
idtaxa_assemblies<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt")
idtaxa_assemblies$type<-"assembly"
idtaxa_reads<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_reads.txt")
idtaxa_reads$type<-"reads"
idtaxa<-rbind(idtaxa_assemblies,idtaxa_reads)

#get family_level assignment

l_tmp<-str_split(idtaxa$assignment,";") 
bac_family<-plyr::ldply(l_tmp, rbind)[7]
idtaxa$label<-bac_family[,1]

idtaxa_occ_target<-idtaxa %>% select(-assignment,-file) %>% 
  filter(bac_family %in% groups_to_include) %>%
  group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))



## 4. get together metaxa and idtaxa
rdna<-rbind(idtaxa_occ_target,metaxa_occ_target)  %>%
pivot_wider(-n,values_from = occurrence, values_fill = 0, names_from = label) %>%
  filter(metagenome %in% mag_presence$metagenome)

write.table(rdna,"analysis/03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv",sep="\t",quote = F, row.names = F)


# add 0 lines for metagenomes that didn't have any of the rDNA 
#only one: X13
assembly_presence<-rdna %>% filter(type=="assembly")
mtg_to_show$metagenome[!(mtg_to_show$metagenome %in% assembly_presence$metagenome)]
X13<-data.frame("type"="assembly","metagenome"="X13","Acetobacteraceae"=0,
"Acidobacteriaceae"=0, "Beijerinckiaceae"=0,  "Sphingomonadaceae"=0,
"Trebouxiophyceae"=0,  "Cystobasidiomycetes"=0, "Tremellomycetes"=0 )
assembly_presence<-rbind(assembly_presence,X13)
  
reads_presence<-rdna %>% filter(type=="reads")





# 3. Visualize
column_order_vector<-c("Acetobacteraceae"=1,"Acidobacteriaceae"=2,
  "Beijerinckiaceae"=3,
  "Sphingomonadaceae"=4,"Trebouxiophyceae"=5,"Cystobasidiomycetes"=6,"Tremellomycetes"=7)
# 3a. mags
M = as.matrix(mag_presence[,2:ncol(mag_presence)])
rownames(M) = mag_presence$metagenome
M  =M[,colSums(M) >0 ]

o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=30, repr.plot.height=10)
node_colors = c("Trebouxiophyceae" = "#00CD6C", 
                "Cystobasidiomycetes" = "#09b1db",                  
                "Tremellomycetes" = "#FFC61E",
                "Acetobacteraceae" = "#c43b0e",
                "Beijerinckiaceae" = "#783f04",
                "Acidobacteriaceae" = "#6c44da",
                "Sphingomonadaceae" = "#34b7eb"
)

#ta = HeatmapAnnotation(df = data.frame(class = clade), col = list(class = node_colors))
cols<-colorRampPalette(brewer.pal(9, "Reds"))(9)
HM = Heatmap(M, show_row_names = F, show_column_names = T, name = " ",
             #top_annotation = ta, 
             column_order = column_order_vector,
             #column_split = clade, 
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
             col = c("0" = "#ffffff", "1" = "#000000")
)
HM



# 3b. assemblies
M1 = as.matrix(assembly_presence[,3:ncol(assembly_presence)])
rownames(M1) = assembly_presence$metagenome
M1  =M1[,colSums(M1) >0 ]

o11 = seriate(dist(M1), method = "TSP")
o12 = seriate(dist(t(M1)), method = "TSP")


HM_assembly = Heatmap(M1, show_row_names = F, show_column_names = T, name = " ",
             #top_annotation = ta, 
             column_order = column_order_vector,
             #column_split = clade, 
             #cluster_rows = row_clust,
             rect_gp = gpar(col= "#cecece"), 
             cluster_rows = reorder(row_clust, get_order(o11)) ,
             row_dend_width = unit(8, "cm"),
             #cell_fun = function(j, i, x, y, width, height, fill) {
             #          grid.rect(x = x, y = y, width = width, height = height, 
             #            gp = gpar(col = 'white', fill = 'white'))
             #          grid.circle(x = x, y = y, r = 0.005,
             #                gp = gpar(fill = fill, col = "#cecece"))
             #    }    ,
             col = c("0" = "#ffffff", "1" = "#000000")
)

HM_assembly


# 3c. reads
M2 = as.matrix(reads_presence[,3:ncol(reads_presence)])
rownames(M2) = reads_presence$metagenome
M2  =M2[,colSums(M2) >0 ]

o21 = seriate(dist(M2), method = "TSP")
o22 = seriate(dist(t(M2)), method = "TSP")


HM_reads = Heatmap(M2, show_row_names = F, show_column_names = T, name = " ",
                      #top_annotation = ta, 
                      column_order = column_order_vector,
                      #column_split = clade, 
                      #cluster_rows = row_clust,
                      rect_gp = gpar(col= "#cecece"), 
                      cluster_rows = reorder(row_clust, get_order(o21)) ,
                      row_dend_width = unit(8, "cm"),
                      #cell_fun = function(j, i, x, y, width, height, fill) {
                      #          grid.rect(x = x, y = y, width = width, height = height, 
                      #            gp = gpar(col = 'white', fill = 'white'))
                      #          grid.circle(x = x, y = y, r = 0.005,
                      #                gp = gpar(fill = fill, col = "#cecece"))
                      #    }    ,
                   col = c("0" = "#ffffff", "1" = "#000000")
)

HM_reads

HM+HM_assembly+HM_reads

pdf(file="analysis/05_MAGs/exploratory_fig/multilevel_screening.pdf",width=10,height=7)
HM+HM_assembly+HM_reads
dev.off()

### 4. make a separate plot for Lichenihabitans
mags_role_filtered$label2<-mags_role_filtered$label
mags_role_filtered$label2[mags_role_filtered$bac_genus=="Lichenihabitans"]<-"Lichenihabitans"
mag_lichenihab<-mags_role_filtered %>% select(metagenome,breadth,label2) %>%
  filter(label2=="Lichenihabitans" | label2=="Mycobiont") %>% 
  filter(metagenome %in% mtg_to_show$metagenome) %>%
  filter(breadth>50) %>% select(-breadth) %>% distinct() %>%
  group_by(metagenome,label2) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label2) %>%
  select(-Mycobiont) %>% mutate(type="MAG",occ=Lichenihabitans) %>% select(-Lichenihabitans)


idtaxa_lichenihab<-idtaxa %>% filter(grepl("Lichenihabitans", assignment)) %>%
     group_by(metagenome,type)  %>% summarize(n=n()) %>% 
  mutate(occ=ifelse(n>0,1,0)) %>% select(-n) %>% filter(metagenome %in% mtg_to_show$metagenome)

licheinhab<-rbind(idtaxa_lichenihab,mag_lichenihab) %>%
  pivot_wider(values_from = occ, names_from = type,values_fill=0)


column_order_vector<-c("MAG"=1,"assembly"=2,"reads"=3)
M3 = as.matrix(licheinhab[,2:ncol(licheinhab)])
rownames(M3) = licheinhab$metagenome


HM_lichenihab = Heatmap(M3, show_row_names = F, show_column_names = T, name = " ",
                   #top_annotation = ta, 
                   column_order = c("MAG","assembly","reads"),
                   #column_split = clade, 
                   #cluster_rows = row_clust,
                   rect_gp = gpar(col= "#cecece"), 
                   cluster_rows = reorder(row_clust, get_order(o21)) ,
                   row_dend_width = unit(8, "cm"),
                   #cell_fun = function(j, i, x, y, width, height, fill) {
                   #          grid.rect(x = x, y = y, width = width, height = height, 
                   #            gp = gpar(col = 'white', fill = 'white'))
                   #          grid.circle(x = x, y = y, r = 0.005,
                   #                gp = gpar(fill = fill, col = "#cecece"))
                   #    }    ,
                   col = c("0" = "#ffffff", "1" = "#000000"))

pdf(file="analysis/05_MAGs/exploratory_fig/multilevel_screening_lichenihabitans.pdf",width=5,height=7)
HM_lichenihab
dev.off()




#only one: X13
assembly_presence<-rdna %>% filter(type=="assembly")
mtg_to_show$metagenome[!(mtg_to_show$metagenome %in% assembly_presence$metagenome)]
X13<-data.frame("type"="assembly","metagenome"="X13","Acetobacteraceae"=0,
                "Acidobacteriaceae"=0, "Beijerinckiaceae"=0,  "Sphingomonadaceae"=0,
                "Trebouxiophyceae"=0,  "Cystobasidiomycetes"=0, "Tremellomycetes"=0 )
assembly_presence<-rbind(assembly_presence,X13)

reads_presence<-rdna %>% filter(type=="reads")


# 4. calculate stats
n_row<-assembly_presence %>% nrow()

rdna_occ_stat<-rdna %>% pivot_longer(-c(type,metagenome),values_to="presence", names_to="lineage") %>%
  group_by(lineage,type) %>% summarize(occ_total=sum(presence)) %>%
  mutate(occ_perc=occ_total*100/n_row) %>% select(-occ_total) %>%
  pivot_wider(names_from=type,values_from=occ_perc)

mag_occ_stat<-mag_presence  %>% pivot_longer(-c(metagenome),values_to="presence", names_to="lineage") %>%
  group_by(lineage) %>% summarize(occ_total=sum(presence)) %>%
  mutate(MAGs=occ_total*100/n_row) %>% select(-occ_total) 

occ_stats<-left_join(rdna_occ_stat,mag_occ_stat)
write.table(occ_stats,"analysis/03_metagenome_reanalysis/occurrence_rDNA_stats.tsv",sep="\t",quote = F, row.names = F)



### for lichenihabitans only
mag_lichenihab %>% filter(occ==1)%>% nrow() / n_row
licheinhab %>% filter(assembly==1)%>% nrow() / n_row
licheinhab %>% filter(reads==1)%>% nrow() / n_row


##make a list of top bacterial families in assemblies and reads. ranked by the 3 of metagenomes they occurred in
#assemblies
top_bacteria_assembly<-idtaxa %>% filter(type=="assembly",!is.na(label)) %>% select(metagenome,label) %>%
  distinct() %>%  group_by(label) %>% summarize(n=n()) %>%  arrange(desc(n))
write.table(top_bacteria_assembly,"analysis/05_MAGs/tables/bacteria_dominant_groups/idtaxa_assembly_bac_families.tsv",sep="\t",quote = F, row.names = F)

#reads
top_bacteria_reads<-idtaxa %>% filter(type=="reads",!is.na(label)) %>% select(metagenome,label) %>%
  distinct() %>%  group_by(label) %>% summarize(n=n()) %>%  arrange(desc(n))
write.table(top_bacteria_reads,"analysis/05_MAGs/tables/bacteria_dominant_groups/idtaxa_reads_bac_families.tsv",sep="\t",quote = F, row.names = F)

