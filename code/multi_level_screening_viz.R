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

##define the list of metagenomes that yeilded at least one mag
metagenomes_with_mags<-mags_role$metagenome %>% unique


##make a table with mag occurrences of the groups of interest
groups_to_include<-c("Trebouxiophyceae","Acidobacteriaceae",
                     "Beijerinckiaceae","Acetobacteraceae",
                     "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes")

mag_presence <- mags_role %>% dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% groups_to_include) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label)

##add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence <- data_frame("metagenome"=metagenomes_with_mags) %>% left_join(mag_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0))

# 2. get info from metaxa for 18S of CYpho, Tremella, and Trebouxia. screening of assemblies and reads
matrix_level5<-read.delim("analysis/03_metagenome_reanalysis/metaxa_level_5_combined.txt")
matrix_level5_long<-matrix_level5 %>% pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% mutate(occurrence=ifelse(abundance>0,1,0))

##summarize by group, reduce the table only to target groups 
metaxa_occ_target<-matrix_level5_long %>% dplyr::filter(occurrence>0, grepl( "Cystobasidiomycetes", Taxa) |
                                                          grepl( "Tremellomycetes", Taxa) |
                                                          grepl( "Trebouxi", Taxa)) %>%
  mutate(label=ifelse(grepl("Cystobasidiomycetes",Taxa),"Cystobasidiomycetes",
                      ifelse(grepl("Beijerinckiaceae",Taxa), "Beijerinckiaceae",
                             ifelse(grepl("Tremellomycetes", Taxa),"Tremellomycetes",
                                    ifelse(grepl("Trebouxi", Taxa),"Trebouxiophyceae",NA))))) %>%
  dplyr::select(-Taxa,-abundance) %>% group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))

# 3. get info from idtaxa for 16S screening
idtaxa_assemblies<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt")
idtaxa_assemblies$type<-"assembly"
idtaxa_reads<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_reads.txt")
idtaxa_reads$type<-"reads"
idtaxa<-rbind(idtaxa_assemblies,idtaxa_reads) %>% filter(metagenome %in% metagenomes_with_mags) #added filtering step to make sure that only intended metagenomes are counted, i.e. the metagenomes with a mycobiont mag and those that were not excluded as misidentified

#get family_level assignment

l_tmp<-str_split(idtaxa$assignment,";") 
bac_family<-plyr::ldply(l_tmp, rbind)[7]
bac_order<- plyr::ldply(l_tmp, rbind)[6]
bac_phylum<-plyr::ldply(l_tmp, rbind)[4]
idtaxa$label<-bac_family[,1]
idtaxa$bac_order<-bac_order[,1]
idtaxa$bac_phylum<-bac_phylum[,1]

idtaxa_occ_target<-idtaxa %>% dplyr::select(-assignment,-file,-bac_order,-bac_phylum) %>% 
  dplyr::filter(label %in% groups_to_include) %>%
  group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))



# 4. get together metaxa and idtaxa
rdna<-rbind(idtaxa_occ_target,metaxa_occ_target)  %>%
  pivot_wider(-n,values_from = occurrence, values_fill = 0, names_from = label) 

write.table(rdna,"analysis/03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv",sep="\t",quote = F, row.names = F)



# 5. Prep. fungal tree for visualization
###load euk tree
tree <- read.tree("analysis/05_MAGs/trees/eukaryotes/Fungi/phylogeny/Concatenated_IQTREE/concat_putative_renamed.contree")


###drop the tips that are not mycobionts. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
mycobiont_mags<-mags_role_filtered %>% dplyr::filter(confirmed_role=="mycobiont") %>% dplyr::select(Genome) %>%
  distinct()
ind <- which(outer( mycobiont_mags$Genome, tree$tip.label, Vectorize(grepl)), arr.ind = T)

mycobiont_tree <- ape::drop.tip(tree, tree$tip.label[-ind[,2]])
#root the tree
mycobiont_tree <- root(mycobiont_tree,outgroup=c("public_SRR14721957_metabat2_bin.1_mycobiont_Ephebe_solida",
                                                 "public_SRR1531569_concoct_bin.26_mycobiont_Peltula_cylindrica",
                                                 "public_SRR14722086_metabat2_bin.14_mycobiont_fungal_sp._Tripp_6058_COLO-L-0051364"),
                       resolve.root = TRUE)



###prepped tree for the heatmap
ta = table(mycobiont_mags)
mtg_to_show<-mags_role_filtered %>% dplyr::filter(confirmed_role=="mycobiont")

modded_tree = mycobiont_tree
###changed tip names to match the Genome column
n1<-str_split_fixed(modded_tree$tip.label, "_", 5) %>% data.frame() %>% mutate(name=paste(X1,X2,X3,X4,sep="_"))
modded_tree$tip.label <- n1$name

for (species in names(ta)){
  leafes_needed  = sum(mtg_to_show$Genome == species)
  while(sum(modded_tree$tip.label == species) < leafes_needed){
    modded_tree = suppressWarnings(bind.tip(modded_tree, species, edge.length=0.0001, where=which(modded_tree$tip.label==species)[1], position=0))
  }
  modded_tree$tip.label[modded_tree$tip.label == species] = mtg_to_show$metagenome[mtg_to_show$Genome == species]
  
}
ape::write.tree(modded_tree,"analysis/05_MAGs/tmptree.tre",digits = 5)
row_clust = ReadDendrogram("analysis/05_MAGs/tmptree.tre",keepRoot = T)


## 6. Visualize
###make a list of metagenomes to visualize: only those with mycobiont mags
metagenomes_with_myco<-mags_role %>% filter(confirmed_role=="mycobiont") %>%
  select(metagenome) %>% distinct()

column_order_vector<-c("Acetobacteraceae"=1,"Acidobacteriaceae"=2,
                       "Beijerinckiaceae"=3,
                       "Sphingomonadaceae"=4,"Trebouxiophyceae"=5,"Cystobasidiomycetes"=6,"Tremellomycetes"=7)
### 6a. mags
mag_presence2<-mag_presence %>% filter(metagenome %in% metagenomes_with_myco$metagenome)
mag_presence2<-mag_presence2[match(rev(labels(row_clust)), mag_presence2$metagenome), ] #reorder

M = as.matrix(mag_presence2[,2:ncol(mag_presence2)])
rownames(M) = mag_presence2$metagenome
M  =M[,colSums(M) >0 ]

HM = Heatmap(M, show_row_names = F, show_column_names = T, name = " ",
             column_order = column_order_vector,
             row_order = rownames(M),
             rect_gp = gpar(col= "#cecece"), 
             row_dend_width = unit(8, "cm"),
             col = c("0" = "#ffffff", "1" = "#000000")
)
HM



# 6b. assemblies
assembly_presence<-rdna %>% dplyr::filter(type=="assembly", metagenome %in% metagenomes_with_myco$metagenome) 
##add empty lines for the metagenomes that didn't have any mags of interest
assembly_presence <- data_frame("metagenome"=metagenomes_with_mags) %>% left_join(assembly_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0))
assembly_presence<-assembly_presence[match(rev(labels(row_clust)), assembly_presence$metagenome), ] #reorder

M1 = as.matrix(assembly_presence[,3:ncol(assembly_presence)])
rownames(M1) = assembly_presence$metagenome
M1  =M1[,colSums(M1) >0 ]

HM_assembly = Heatmap(M1, show_row_names = F, show_column_names = T, name = " ",
                      column_order = column_order_vector,
                      row_order = rownames(M1),
                      rect_gp = gpar(col= "#cecece"), 
                      row_dend_width = unit(8, "cm"),
                      col = c("0" = "#ffffff", "1" = "#000000")
)

HM_assembly


# 6c. reads
reads_presence<-rdna %>% dplyr::filter(type=="reads", metagenome %in% metagenomes_with_myco$metagenome) 
##add empty lines for the metagenomes that didn't have any mags of interest
reads_presence <- data_frame("metagenome"=metagenomes_with_mags) %>% left_join(reads_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0))

reads_presence<-reads_presence[match(rev(labels(row_clust)), reads_presence$metagenome), ] #reorder

M2 = as.matrix(reads_presence[,3:ncol(reads_presence)])
rownames(M2) = reads_presence$metagenome


###add right annotation showing mycobiont taxonomy
mtg_info <- read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
mtg_info$label<-mtg_info$order
mtg_info$label[is.na(mtg_info$label)]<-"Unknown"

lichen_group<-mtg_info %>% dplyr::select(Run,label)
mycobiont_group_df<-data.frame(rownames(M2)) %>% left_join(lichen_group,by=c("rownames.M2."="Run"))
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

ra = rowAnnotation(df=data.frame(mycobiont_group = mycobiont_group),
                   col = list(mycobiont_group = myco_node_color),
                   annotation_legend_param = list(mycobiont_group=list(title = "LFS order")))


HM_reads = Heatmap(M2, show_row_names = F, show_column_names = T, name = " ",
                   column_order = column_order_vector,
                   row_order = rownames(M2),
                   right_annotation = ra, 
                   rect_gp = gpar(col= "#cecece"), 
                   row_dend_width = unit(8, "cm"),
                   col = c("0" = "#ffffff", "1" = "#000000")
)

HM_reads


svglite(file="analysis/05_MAGs/exploratory_fig/multilevel_screening.svg",width=6,height=8)
HM+HM_assembly+HM_reads
dev.off()


### 7. make a separate plot for Lichenihabitans
mags_role_filtered$label2<-mags_role_filtered$label
mags_role_filtered$label2[mags_role_filtered$bac_genus=="Lichenihabitans"]<-"Lichenihabitans"
mag_lichenihab<-mags_role_filtered %>% dplyr::select(metagenome,breadth,label2) %>%
  dplyr::filter(label2=="Lichenihabitans") %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,label2) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label2) %>%
 mutate(occ=Lichenihabitans) %>% dplyr::select(-Lichenihabitans)
##add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_lichenihab <- data_frame("metagenome"=metagenomes_with_mags) %>% left_join(mag_lichenihab) %>%
  mutate_if(is.numeric, ~replace_na(.,0)) %>%  mutate(type="MAG")

idtaxa_lichenihab<-idtaxa %>% dplyr::filter(grepl("Lichenihabitans", assignment)) %>%
     group_by(metagenome,type)  %>% summarize(n=n()) %>% 
  mutate(occ=ifelse(n>0,1,0)) %>% dplyr::select(-n) %>% dplyr::filter(metagenome %in% metagenomes_with_mags)

licheinhab<-rbind(idtaxa_lichenihab,mag_lichenihab) %>%
  pivot_wider(values_from = occ, names_from = type,values_fill=0)


column_order_vector<-c("MAG"=1,"assembly"=2,"reads"=3)
licheinhab2<-licheinhab %>% filter(metagenome %in% metagenomes_with_myco$metagenome)
M3 = as.matrix(licheinhab2[,2:ncol(licheinhab2)])
rownames(M3) = licheinhab2$metagenome


HM_lichenihab = Heatmap(M3, show_row_names = F, show_column_names = T, name = " ",
                   #top_annotation = ta, 
                   column_order = c("MAG","assembly","reads"),
                   #column_split = clade, 
                   #cluster_rows = row_clust,
                   rect_gp = gpar(col= "#cecece"), 
                   cluster_rows = rev(row_clust),
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





# 8. calculate stats
n_mtg<-metagenomes_with_mags %>% length()

rdna_occ_stat<-rdna %>% filter(metagenome %in% metagenomes_with_mags) %>%
  pivot_longer(-c(type,metagenome),values_to="presence", names_to="lineage") %>%
  group_by(lineage,type) %>% summarize(occ_total=sum(presence)) %>%
  mutate(occ_perc=occ_total*100/n_mtg) %>% dplyr::select(-occ_total) %>%
  pivot_wider(names_from=type,values_from=occ_perc)

mag_occ_stat<-mag_presence  %>% pivot_longer(-c(metagenome),values_to="presence", names_to="lineage") %>%
  group_by(lineage) %>% summarize(occ_total=sum(presence)) %>%
  mutate(MAGs=occ_total*100/n_mtg) %>% dplyr::select(-occ_total) 

occ_stats<-left_join(rdna_occ_stat,mag_occ_stat)
write.table(occ_stats,"analysis/03_metagenome_reanalysis/occurrence_rDNA_stats.tsv",sep="\t",quote = F, row.names = F)



### for lichenihabitans only
mag_lichenihab %>% dplyr::filter(occ==1)%>% nrow() / n_mtg
licheinhab %>% dplyr::filter(assembly==1)%>% nrow() / n_mtg
licheinhab %>% dplyr::filter(reads==1)%>% nrow() / n_mtg


##make a list of top bacterial families in assemblies and reads. ranked by the % of metagenomes they occurred in
#assemblies
top_bacteria_assembly<-idtaxa %>% dplyr::filter(type=="assembly",!is.na(label),metagenome %in% metagenomes_with_mags) %>% 
  dplyr::select(metagenome,label) %>%
  dplyr::distinct() %>%  dplyr::group_by(label) %>% dplyr::summarize(n=n()) %>%
  mutate(percentage=n*100/n_mtg) %>%
dplyr::arrange(dplyr::desc(n))
write.table(top_bacteria_assembly,"analysis/05_MAGs/tables/bacteria_dominant_groups/idtaxa_assembly_bac_families.tsv",sep="\t",quote = F, row.names = F)

#reads
top_bacteria_reads<-idtaxa %>% dplyr::filter(type=="reads",!is.na(label)) %>% 
  dplyr::select(metagenome,label) %>%
  dplyr::distinct() %>%  dplyr::group_by(label) %>%
  dplyr::summarize(n=n()) %>%  
  mutate(percentage=n*100/n_mtg) %>%
  dplyr::arrange(dplyr::desc(n))
write.table(top_bacteria_reads,"analysis/05_MAGs/tables/bacteria_dominant_groups/idtaxa_reads_bac_families.tsv",sep="\t",quote = F, row.names = F)

#combine both tables to make a supplementary table
top_idtaxa<-rbind(top_bacteria_assembly %>% mutate(Source="Assemblies"),top_bacteria_reads %>% mutate(Source="Reads")) %>%
          mutate(Number_of_metagenomes_detected = n, Prevalence = n*100/(mag_presence$metagenome %>% unique %>% length)) %>%
            left_join(idtaxa %>% select(label,bac_order,bac_phylum) %>% distinct()) %>%
          select(Source,label,bac_order,bac_phylum,Prevalence,Number_of_metagenomes_detected)
write.table(top_idtaxa,"results/tables/idtaxa_top_bac_families.tsv",sep="\t",quote = F, row.names = F)



