#setwd("~/Documents/phd/coverage")
library(tidyverse)
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
source("code/utils.R")

#select high quality euk mags for annotation

#load table
mags<-read.delim("analysis/05_MAGs/tables/MAG_suppl_table.tsv")

#filtered and got 84 mags
mags_sel<-mags %>% filter(completeness>=95,contamination<10,!grepl("d__Bacteria",classification))

mags_sel<-mags_sel %>% mutate(group=ifelse(grepl("chloro",classification),"alga","fungus")) %>%
  mutate(locustag = paste0(str_match(Genome,"^[a-z]*_([a-zA-Z0-9]*)_.*")[,2],"b",str_match(Genome,"^.*\\.([0-9]*)$")[,2])) %>%
  select(Genome,group,locustag)

#write
write.table(mags_sel,"analysis/12_annotate_euks/euk_annot_list.txt",sep="\t",quote = F, row.names = F,col.names = T)

#add info
tax_info<-read.delim("submission/mag_submission/euk_ncbi.txt")
mags_sel2 <- mags_sel %>% left_join(tax_info) %>% select(-classification) %>%
  dplyr::rename( classification = classification_fixed, 
          role_in_lichen=confirmed_role,
         lichen_species=Lichen.metagenomes)
write.table(mags_sel2,"analysis/12_annotate_euks/euk_annot_table.txt",sep="\t",quote = F, row.names = F,col.names = T)

  

#which metagenomes had high quality mycobiont mag + high quality photobiont mag + at least one high quality bacterial mag?
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
complete <- mags_role %>% left_join(mags) %>% filter(completeness>=95) %>%
  group_by(metagenome,confirmed_role) %>% summarize(n=n()) %>%
  pivot_wider(names_from = confirmed_role, values_from = n,values_fill = 0) %>%
  filter(mycobiont>0 & bacteria_other>0 & (photobiont_chloro>0 | photobiont_cyano>0))
#in total had 22 metagenomes there

#which metagenomes had high quality mycobiont mag + high quality photobiont mag + at least one of the 63 annotated mags
prok_ann<-read.delim("analysis/07_annotate_MAGs/mag_table.txt")
complete2 <- mags_role %>% left_join(mags) %>% filter(completeness>=95) %>%
  mutate(confirmed_role2 = ifelse(Genome %in% prok_ann$mag,"bacteria_annot",confirmed_role)) %>%
  group_by(metagenome,confirmed_role2) %>% summarize(n=n()) %>%
  pivot_wider(names_from = confirmed_role2, values_from = n,values_fill = 0) %>%
  filter(mycobiont>0 & bacteria_annot>0 & (photobiont_chloro>0 | photobiont_cyano>0))
#in total had 17 metagenomes there, half from Philpp's paper

#save as a table
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
complete_info<-mtg_info %>% filter(Run %in% complete$metagenome) %>% 
  mutate(contains_annot_bacterial_mags = ifelse(Run %in% complete2$metagenome, T,F)) %>%
  select(Run,Source,Lichen.metagenomes,family,order,class,contains_annot_bacterial_mags)
  
write.table(complete_info,"analysis/12_annotate_euks/mtg_high_qual_mags.txt",sep="\t",quote = F, row.names = F,col.names = T)

#which genomes are present in the 17 metagenomes?
genomes_in_17<-mags_role %>% left_join(mags) %>% filter(completeness>=95,metagenome %in% complete2$metagenome)

#KEGG results
module_info1<-read.delim("analysis/12_annotate_euks/module_descriptions_prokaryotes.txt")
module_info2<-read.delim("analysis/12_annotate_euks/module_descriptions_eukaryotes.txt")
module_info<-rbind(module_info2,module_info1) %>% distinct()

prok_kegg<-read.delim("analysis/12_annotate_euks/hdf_prokaryotes.txt")
euk_kegg<-read.delim("analysis/12_annotate_euks/hdf_eukaryotes.txt")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")

prok_kegg_long<-prok_kegg  %>% pivot_longer(-X, names_to = "Genome",values_to = "completeness") 
euk_kegg_long<-euk_kegg %>%  pivot_longer(-X, names_to = "Genome",values_to = "completeness")
kegg <- prok_kegg_long %>% rbind(euk_kegg_long) %>% left_join(module_info) %>% 
  mutate(module=gsub("\\_.*","",X))

#select biotin and thiamine pathways and only genomes from 17 complete mtg. only keep one definition of module per genome (select the highest completeness)
kegg_sel <- kegg %>% filter(Genome %in% genomes_in_17$Genome, 
            module %in% c("M00899","M00898","M00897",
                     "M00123","M00950","M00573","M00577")) %>%
  group_by(Genome,module) %>% mutate(mean_completeness=mean(completeness)) %>%
  mutate(vitamin=if_else(module %in% c("M00899","M00898","M00897"),"thiamine","biotin"),
       completeness2 = if_else(mean_completeness==1,"Complete","50%-99% complete")) %>%
  select(-c(X,completeness)) %>% distinct()
  
kegg_sel$Genome<-as.factor(kegg_sel$Genome)
kegg_sel$module<-as.factor(kegg_sel$module)

#join with the mags_role table to assign MAG roles
kegg_sel<-kegg_sel %>% inner_join(mags_role %>% select(Genome,bac_family,metagenome,confirmed_role) %>% 
        filter(metagenome %in% complete2$metagenome),relationship = "many-to-many") %>%
  mutate(role = case_when(confirmed_role=="mycobiont" ~"LFS",
                          confirmed_role== "bacteria_other" ~ "Bacteria",
                          confirmed_role=="photobiont_cyano" ~ "Cyanobacterial\nPhotobiont",
                          confirmed_role=="cephalodia_cyano" ~ "Cyanobacterial\nPhotobiont",
                          confirmed_role=="photobiont_chloro" ~ "Alga",
                          confirmed_role=="fungi_other" ~ "Other Fungi")) %>%
  mutate(taxonomy = ifelse(is.na(bac_family),role,bac_family))


kegg_sel$role<-factor(kegg_sel$role,levels=c("LFS","Alga","Cyanobacterial\nPhotobiont","Bacteria","Other Fungi"))
cols2<-c(cols,"LFS"="#522f01","Alga"="#096103","Cyanobacterial\nPhotobiont"="#02819e","Other Fungi"="#180061")


plotting<-function(x){
  ggplot(kegg_sel %>% filter(metagenome== x))+
         geom_point(aes(x=Genome,y=module,shape=completeness2,alpha = I(ifelse(mean_completeness<0.5, 0, 1)),color=taxonomy),size=2)+
  facet_grid(vitamin~role,scales = "free")+
  scale_shape_manual(values=c(16,10),limits = c("Complete","50%-99% complete"))+
  scale_color_manual(values=cols2)+ggtitle(x)+guides(color="none")+
  theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        plot.title = element_text(size=8),
        axis.text.y = element_text(size=6),
        strip.text.x = element_text(size =7),
        strip.text.y = element_text(size =7))}

plotting_no_axis_label<-function(x){
  ggplot(kegg_sel %>% filter(metagenome== x))+
    geom_point(aes(x=Genome,y=module,shape=completeness2,alpha = I(ifelse(mean_completeness<0.5, 0, 1)),color=taxonomy),size=2)+
    facet_grid(vitamin~role,scales = "free")+
    scale_shape_manual(values=c(16,10),limits = c("Complete","50%-99% complete"))+
    scale_color_manual(values=cols2)+ggtitle(x)+guides(color="none")+
    theme_bw()+
    theme(axis.title = element_blank(),axis.text = element_blank(),
          plot.title = element_text(size=8),
          strip.text.x = element_text(size =7),
          strip.text.y = element_text(size =7))}



p1 <- plotting_no_axis_label("X3")
p2 <- plotting("X8")
p3 <- plotting_no_axis_label("X14")
p4 <- plotting("VT1")
p5 <- plotting_no_axis_label("VT12")
p6 <- plotting("SRR5808930")
p7 <- plotting_no_axis_label("T1914")
p8 <- plotting_no_axis_label("T1888")
p9 <- plotting_no_axis_label("SRR11456913")
p10 <- plotting_no_axis_label("SRR11456919")
p11 <- plotting("ERR4179390")
p12 <- plotting_no_axis_label("T1894")
p13 <- plotting_no_axis_label("T1868")
p14 <- plotting("T1916")
p15 <- plotting_no_axis_label("T1904")
p16 <- plotting("T1889")
p17 <- plotting("T1867")

#NB: excluded T1916, because it actually doesn't have the true photobiont MAG. Gyalectas supposed to have Trentepolia, and this metagenome only has a Trebouxia MAG

library(patchwork)                      
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+ guide_area()+plot_layout(ncol = 3,guides = 'collect')
ggsave("results/figures/pathway_cooc.pdf",width=300,height=300,device="pdf",bg="white",units="mm")


p <- p6+p8+p7 + p4+p5+p1 + p2+p3+p15 + p11+p9+p10 + p17+p13 + p12+p16

p+guide_area()+plot_layout(guides = 'collect',ncol = 3)
ggsave("results/figures/pathway_cooc2.pdf",width=200,height=225,device="pdf",bg="white",units="mm")



#how many unannotated mags have these metagenomes?
mags_role %>% mutate(annotated = ifelse(Genome %in% kegg$Genome,"true","false")) %>%
  group_by(metagenome,annotated) %>% summarise(n=n()) %>%
  pivot_wider(names_from = annotated,values_from=n,values_fill = 0) %>%
  mutate(total=true+false) %>%
  filter(metagenome %in% complete2$metagenome) 


#figure with phylogeny and pathways
##fungal tree
library(ape)
treef <- read.tree("analysis/05_MAGs/trees/david_trees_euk/Fungi/phylogeny/Concatenated_IQTREE/concat.contree")

#drop the tips that are not annotated. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
f_id <- which(outer( mags_sel2$Genome, treef$tip.label, Vectorize(grepl)), arr.ind = T)
treef2 <- ape::drop.tip(treef, treef$tip.label[-f_id[,2]])
#root the tree
treef2 <- root(treef2,outgroup=c("public_SRR1531569_concoct_bin.26",
                                                 "public_SRR14722086_metabat2_bin.14"),
                       resolve.root = TRUE)
treef2<-ladderize(treef2,right=F)

#got node order
is_tip <- treef2$edge[,2] <= length(treef2$tip.label)
ordered_tips <- treef2$edge[is_tip, 2]
tip_order_f<-treef2$tip.label[ordered_tips]

write.tree(treef2, "analysis/12_annotate_euks/fungal.tree")
row_clust <- ReadDendrogram("analysis/12_annotate_euks/fungal.tree")


#prep matrix
df_f <- kegg %>% filter(Genome %in% treef2$tip.label, 
                        module %in% c("M00899","M00898","M00897","M00123","M00950","M00573","M00577","M00925","M00122")) %>%
  group_by(Genome,module) %>% mutate(mean_completeness=mean(completeness)) %>%
  mutate(vitamin=case_when(module %in% c("M00899","M00898","M00897")~"thiamine",
                           module %in% c("M00123","M00950","M00573","M00577")~"biotin",
                           T~"cobalamin")) %>%
  select(-c(X,completeness)) %>% distinct()

heatmapf <- df_f %>% select(Genome, module, mean_completeness) %>% 
  pivot_wider(names_from = module,values_from = mean_completeness) %>%
  arrange(match(Genome, tip_order_f))

rownames_heatmapf <- heatmapf$Genome
heatmapf <- as.matrix(heatmapf[,2:ncol(heatmapf)])
rownames(heatmapf) <- rownames_heatmapf 
o1 <- seriate(dist(heatmapf), method = "TSP")

# color function
col_fun = circlize::colorRamp2(c(0, 1),
                               c("white", scales::muted("blue")))

# cluster rows by vitamin
cluster_col <- c("thiamine","thiamine","thiamine","biotin","biotin","biotin","biotin","cobalamin","cobalamin")
names(cluster_col) = colnames(heatmapf)

ht1 <- Heatmap(heatmapf,
               show_column_names = TRUE,
               show_row_names = TRUE,
               col=col_fun,
               column_split=cluster_col,
               cluster_columns = FALSE,
               #column_labels = column_labels[colnames(heatmap)],
               column_names_gp = grid::gpar(fontsize = 4),
               row_names_gp = grid::gpar(fontsize = 8),
               #cluster_rows = reorder(row_clust, get_order(o1)),
               row_order = rownames(heatmapf),
               cluster_rows = F,
               row_title = NULL,
               column_title = NULL,
               width = unit(7, "cm"),
               height = unit(17, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE
)

pdf(file="results/module_completeness_fungi.pdf",width=3, height=6.7)
draw(ht1, heatmap_legend_side = "left")
dev.off()

#get a version of heatmap with the tree to merge with the version below (has to do that because ComplexHeatmaps is buggy)

ht1_2 <- Heatmap(heatmapf,
               show_column_names = TRUE,
               show_row_names = TRUE,
               col=col_fun,
               column_split=cluster_col,
               cluster_columns = FALSE,
               #column_labels = column_labels[colnames(heatmap)],
               column_names_gp = grid::gpar(fontsize = 4),
               row_names_gp = grid::gpar(fontsize = 8),
               cluster_rows = reorder(row_clust, get_order(o1)),
               row_title = NULL,
               column_title = NULL,
               width = unit(7, "cm"),
               height = unit(17, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE
)
pdf(file="results/module_completeness_fungi_tree.pdf",width=1.57, height=6.7)
draw(ht1_2, heatmap_legend_side = "left")
dev.off()

# algal heatmap
treea <- read.tree("analysis/05_MAGs/trees/david_trees_euk/Algae/phylogeny/Concatenated_IQTREE/concat.contree")

#drop the tips that are not annotated. had to use regular expressions, because the new tip labels include metagdata and aren't identical to the mag names
a_id <- which(outer( mags_sel2$Genome, treea$tip.label, Vectorize(grepl)), arr.ind = T)
treea2 <- ape::drop.tip(treea, treea$tip.label[-a_id[,2]])
#root the tree
treea2 <- root(treea2,outgroup=c("public_ERR4179390_concoct_bin_45",
                                 "public_ERR4179390_concoct_bin_75"),
               resolve.root = TRUE)
treea2<-ladderize(treea2,right=F)

#fix names
treea2$tip.label<-str_replace(treea2$tip.label,"bin_","bin.")
treea2$tip.label<-str_replace(treea2$tip.label,"merged_","merged.")

#got node order
is_tip <- treea2$edge[,2] <= length(treea2$tip.label)
ordered_tips <- treea2$edge[is_tip, 2]
tip_order_a<-treea2$tip.label[ordered_tips]

write.tree(treea2, "analysis/12_annotate_euks/algal.tree")
row_clust_a <- ReadDendrogram("analysis/12_annotate_euks/algal.tree")


#prep matrix
df_a <- kegg %>% filter(Genome %in% treea2$tip.label, 
                        module %in% c("M00899","M00898","M00897","M00123","M00950","M00573","M00577","M00925","M00122")) %>%
  group_by(Genome,module) %>% mutate(mean_completeness=mean(completeness)) %>%
  mutate(vitamin=case_when(module %in% c("M00899","M00898","M00897")~"thiamine",
                           module %in% c("M00123","M00950","M00573","M00577")~"biotin",
                           T~"cobalamin")) %>%
  select(-c(X,completeness)) %>% distinct()

heatmapa <- df_a %>% select(Genome, module, mean_completeness) %>% 
  pivot_wider(names_from = module,values_from = mean_completeness) %>%
  arrange(match(Genome, tip_order_a))

rownames_heatmapa <- heatmapa$Genome
heatmapa <- as.matrix(heatmapa[,2:ncol(heatmapa)])
rownames(heatmapa) <- rownames_heatmapa
o2 <- seriate(dist(heatmapa), method = "TSP")

ht2 <- Heatmap(heatmapa,
               show_column_names = TRUE,
               show_row_names = TRUE,
               col=col_fun,
               column_split=cluster_col,
               cluster_columns = FALSE,
               #column_labels = column_labels[colnames(heatmap)],
               column_names_gp = grid::gpar(fontsize = 4),
               row_names_gp = grid::gpar(fontsize = 8),
               #cluster_rows = reorder(row_clust, get_order(o1)),
               row_order = rownames(heatmapa),
               cluster_rows = F,
               row_title = NULL,
               column_title = NULL,
               width = unit(7, "cm"),
               height = unit(3, "cm"),
               heatmap_legend_param = list(
                 legend_direction = "vertical", 
                 legend_width = unit(5, "cm")
               ),
               rect_gp = gpar(col = "white", lwd = 0),
               border=TRUE
)

pdf(file="results/module_completeness_alga.pdf",width=1.57, height=1.18)
draw(ht2, heatmap_legend_side = "left")
dev.off()

#get a version of heatmap with the tree to merge with the version below (has to do that because ComplexHeatmaps is buggy)

ht2_2 <- Heatmap(heatmapa,
                 show_column_names = TRUE,
                 show_row_names = TRUE,
                 col=col_fun,
                 column_split=cluster_col,
                 cluster_columns = FALSE,
                 #column_labels = column_labels[colnames(heatmap)],
                 column_names_gp = grid::gpar(fontsize = 4),
                 row_names_gp = grid::gpar(fontsize = 8),
                 cluster_rows = reorder(row_clust_a, get_order(o2)),
                 row_title = NULL,
                 column_title = NULL,
                 width = unit(7, "cm"),
                 height = unit(3, "cm"),
                 heatmap_legend_param = list(
                   legend_direction = "vertical", 
                   legend_width = unit(5, "cm")
                 ),
                 rect_gp = gpar(col = "white", lwd = 0),
                 border=TRUE
)
pdf(file="results/module_completeness_alga_tree.pdf",width=1.57, height=1.18)
draw(ht2_2, heatmap_legend_side = "left")
dev.off()



# Check InterPro annotation to confirm that they are consistent with KEGG
ipr_target<-read.delim("analysis/12_annotate_euks/KEGG_interpro.txt",sep="\t")

#check manually one genome
VT12_ipr<-read.delim("analysis/12_annotate_euks/gulya_ips_annot/ips_euks/private_X14_concoct_bin.2.annotations.txt",sep="\t") %>%
  select(GeneID,Product,InterPro)
#split InterPro annotations
VT12_ipr2 <-VT12_ipr %>% 
  mutate(InterPro = strsplit(InterPro, ";I")) %>%
  unnest(InterPro) %>% mutate(InterPro=str_replace(InterPro,"^PR","IPR"))
VT12_ipr2$InterProID = gsub( " .*$", "", VT12_ipr2$InterPro)

#any annotations from biotin pathway?
VT12_ipr2 %>% filter(InterProID %in% ipr_target$InterPro) %>% nrow()

#apply to all euks
genome_list<- list.files("analysis/12_annotate_euks/gulya_ips_annot/ips_euks/", pattern="*annotations.txt",full.names=F)

screen_ips<-function(x){
  ipr<-read.delim(paste0("analysis/12_annotate_euks/gulya_ips_annot/ips_euks/",x),sep="\t") %>% select(GeneID,Product,InterPro)
  #split InterPro annotations
  ipr2 <-ipr %>% 
    mutate(InterPro = strsplit(InterPro, ";I")) %>%
    unnest(InterPro) %>% mutate(InterPro=str_replace(InterPro,"^PR","IPR"))
  ipr2$InterProID = gsub( " .*$", "", ipr2$InterPro)
  #any annotations from biotin pathway?
  out<-ipr2 %>% filter(InterProID %in% ipr_target$InterPro) %>% nrow()
  return(data.frame("genome"=x,"biotin_genes"=out))
}  

l<-lapply(genome_list, screen_ips)
biotin_table<-do.call(rbind,l)
biotin_table$genome<-str_replace(biotin_table$genome,".annotations.txt","")
#add mag info
biotin_table<-biotin_table %>% 
  left_join(mags_sel2 %>% select(Genome,classification),by=c("genome"="Genome"))

