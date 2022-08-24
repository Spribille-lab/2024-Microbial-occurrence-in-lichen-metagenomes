# 1. Make a table for presence/absence of mags

#load the mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")


#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"


#groups to include in the heatmap
groups_to_include<-c("Acetobacteraceae","Caulobacteraceae",
                     "Ktedonobacteraceae","UBA10450",
                     "Burkholderiaceae","Capsulimonadaceae","Streptosporangiaceae")

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




## 3. get info from idtaxa for 16S screening
idtaxa_assemblies<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt")
idtaxa_assemblies$type<-"assembly"
idtaxa_reads<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_reads.txt")
idtaxa_reads$type<-"reads"
idtaxa<-rbind(idtaxa_assemblies,idtaxa_reads)


#get family_level assignment

l_tmp<-str_split(idtaxa$assignment,";") 
bac_family<-plyr::ldply(l_tmp, rbind)[7]
idtaxa$bac_family<-bac_family[,1]

idtaxa_occ_target<-idtaxa %>% select(-assignment,-file) %>% 
 filter(bac_family %in% groups_to_include) %>%
  group_by(type,metagenome,bac_family) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))


rdna<-idtaxa_occ_target  %>%
  pivot_wider(-n,values_from = occurrence, values_fill = 0, names_from = bac_family) %>%
  filter(metagenome %in% mag_presence$metagenome)



# add 0 lines for metagenomes that didn't have any of the rDNA 
#only one: X13
assembly_presence<-rdna %>% filter(type=="assembly")
empty_lines<-mtg_to_show$metagenome[!(mtg_to_show$metagenome %in% assembly_presence$metagenome)]
empty_lines_df<-data.frame("type"="assembly","metagenome"=empty_lines,"UBA10450"=rep(0,length(empty_lines)),
                "Caulobacteraceae"=rep(0,length(empty_lines)), "Ktedonobacteraceae"=rep(0,length(empty_lines)),  "Burkholderiaceae"=rep(0,length(empty_lines)),
                "Acetobacteraceae"=rep(0,length(empty_lines)),  "Capsulimonadaceae"=rep(0,length(empty_lines)), "Streptosporangiaceae"=rep(0,length(empty_lines)) )
assembly_presence<-rbind(assembly_presence,empty_lines_df)

reads_presence<-rdna %>% filter(type=="reads")
empty_lines2<-mtg_to_show$metagenome[!(mtg_to_show$metagenome %in% reads_presence$metagenome)]
empty_lines2_df<-data.frame("type"="reads","metagenome"=empty_lines2,"UBA10450"=rep(0,length(empty_lines2)),
                           "Caulobacteraceae"=rep(0,length(empty_lines2)), "Ktedonobacteraceae"=rep(0,length(empty_lines2)),  "Burkholderiaceae"=rep(0,length(empty_lines2)),
                           "Acetobacteraceae"=rep(0,length(empty_lines2)),  "Capsulimonadaceae"=rep(0,length(empty_lines2)), "Streptosporangiaceae"=rep(0,length(empty_lines2)) )
reads_presence<-rbind(reads_presence,empty_lines2_df)




# 3. Visualize
column_order_vector<-c("UBA10450"=1,"Caulobacteraceae"=2,
                       "Ktedonobacteraceae"=3,
                       "Burkholderiaceae"=4,"Acetobacteraceae"=5,"Capsulimonadaceae"=6,"Streptosporangiaceae"=7)
# 3a. mags
M = as.matrix(mag_presence[,2:ncol(mag_presence)])
rownames(M) = mag_presence$metagenome
M  =M[,colSums(M) >0 ]

o1 = seriate(dist(M), method = "TSP")
o2 = seriate(dist(t(M)), method = "TSP")

#plot
options(repr.plot.width=30, repr.plot.height=10)
node_colors = c("UBA10450" = "#00CD6C", 
                "Caulobacteraceae" = "#09b1db",                  
                "Ktedonobacteraceae" = "#FFC61E",
                "Burkholderiaceae" = "#c43b0e",
                "Capsulimonadaceae" = "#783f04",
                "Streptosporangiaceae" = "#6c44da",
                "Acetobacteraceae" = "black"
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

pdf(file="analysis/05_MAGs/exploratory_fig/multilevel_screening_other_bacteria.pdf",width=10,height=7)
HM+HM_assembly+HM_reads
dev.off()

