source("code/utils.R")
library(tidyverse)
library(svglite)
library(cowplot)
library(patchwork)


# 1. Make a table for presence/absence of mags. include only metagenomes that yielded at least one mag

##load the mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

##add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Trebouxiophyceae"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cystobasidiomycetes"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremellomycetes"
mags_role$label[mags_role$bac_phylum=="Cyanobacteria"]<-"Cyanobacteria"
mags_role<-mags_role %>%mutate(bac_family2=ifelse(bac_family=="Unknown",paste(bac_order," fam.",sep=""),bac_family))
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family2," gen. sp.",sep=""),bac_genus))
mags_role$bac_genus2[mags_role$bac_genus2=="Acetobacteraceae gen. sp." & !is.na(mags_role$bac_genus2)]<-"unclassified Acetobacteraceae"
#mags_role$label[mags_role$Genome %in% c("public_SRR1532736_concoct_bin.20","private_GTX0163_concoct_bin.46",
 #           "public_SRR13125477_metabat2_bin.13","private_X2_concoct_bin.1",
  #          "private_T1882_metabat2_merged.0","private_X11_concoct_bin.15",
   #         "private_T1912_concoct_bin.0","private_X10_concoct_bin.7",
    #        "public_SRR14722116_metabat2_bin.5","public_SRR14722120_metabat2_bin.4",
     #       "public_SRR7232214_concoct_bin.9","public_SRR7232212_concoct_bin.0",
      #      "public_SRR7232213_concoct_bin.0","private_GTX0161_concoct_bin.0",
       #     "private_TS1974_concoct_bin.29","public_SRR2387885_metabat2_bin.2",
        #    "public_SRR7232211_concoct_bin.8","public_SRR5808932_metabat2_bin.2",
         #   "private_X1_metabat2_bin.11","private_T1904_concoct_bin.38")]<-"Trebouxiaceae"
#mags_role$label[mags_role$Genome %in% c("public_SRR11456915_concoct_bin.5",
 #           "public_SRR13685159_concoct_merged.0","public_SRR11456919_concoct_bin.22",
  #          "public_SRR14722188_concoct_bin.30","public_SRR14722042_metabat2_bin.4",
   #         "private_T1894_concoct_bin.66","public_SRR14722152_metabat2_bin.2",
    #        "private_T1889_concoct_bin.16","public_SRR14722135_metabat2_bin.5")]<-"Coccomyxa-Elliptochloris clade"

##define the list of metagenomes that yeilded a mycobiont MAG
metagenomes_with_myco<-mags_role$metagenome[mags_role$confirmed_role=="mycobiont"] %>% unique


##make a table with mag occurrences of the groups of interest
groups_to_include<-c("Trebouxiophyceae","Acidobacteriaceae",
                    "Beijerinckiaceae","Acetobacteraceae",
                    "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes","Cyanobacteria")
#groups_to_include<-c("Acidobacteriaceae","Trebouxiaceae","Coccomyxa-Elliptochloris clade",
                   # "Beijerinckiaceae","Acetobacteraceae",
                  # "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes","Cyanobacteria")



mag_presence <- mags_role %>% dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% groups_to_include) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label)

##add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence <- data_frame("metagenome"=metagenomes_with_myco) %>% left_join(mag_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0))

# 2. get info from metaxa for 18S of CYpho, Tremella, Ulvophyceae, and Trebouxia. screening of assemblies and reads
matrix_level5<-read.delim("analysis/03_metagenome_reanalysis/metaxa_level_5_combined.txt")
matrix_level5_long<-matrix_level5 %>% pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% mutate(occurrence=ifelse(abundance>0,1,0))

##summarize by group, reduce the table only to target groups 
metaxa_occ_target<-matrix_level5_long %>% dplyr::filter(occurrence>0, grepl( "Cystobasidiomycetes", Taxa) |
                                                          grepl( "Tremellomycetes", Taxa) |
                                                          grepl( "Trebouxi", Taxa) |
                                                          grepl( "Ulvophyc", Taxa)) %>%
  mutate(label=ifelse(grepl("Cystobasidiomycetes",Taxa),"Cystobasidiomycetes",
                      ifelse(grepl("Beijerinckiaceae",Taxa), "Beijerinckiaceae",
                             ifelse(grepl("Tremellomycetes", Taxa),"Tremellomycetes",
                                    ifelse(grepl("Trebouxi", Taxa),"Trebouxiophyceae",
                                           ifelse(grepl( "Ulvophyc", Taxa),"Ulvophyceae",NA)))))) %>%
  dplyr::select(-Taxa,-abundance) %>% group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))

# 3. get info from idtaxa for 16S screening
idtaxa_assemblies<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt")
idtaxa_assemblies$type<-"assembly"
idtaxa_reads<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_reads.txt")
idtaxa_reads$type<-"reads"
idtaxa<-rbind(idtaxa_assemblies,idtaxa_reads) %>% filter(metagenome %in% metagenomes_with_myco) #added filtering step to make sure that only intended metagenomes are counted, i.e. the metagenomes with a mycobiont mag and those that were not excluded as misidentified

#get family_level assignment

l_tmp<-str_split(idtaxa$assignment,";") 
bac_family<-plyr::ldply(l_tmp, rbind)[7]
bac_order<- plyr::ldply(l_tmp, rbind)[6]
bac_phylum<-plyr::ldply(l_tmp, rbind)[4]
idtaxa$label<-bac_family[,1]
idtaxa$bac_order<-bac_order[,1]
idtaxa$bac_phylum<-bac_phylum[,1]
idtaxa$label[idtaxa$bac_phylum=="Cyanobacteria"]<-"Cyanobacteria"

idtaxa_occ_target<-idtaxa %>% dplyr::select(-assignment,-file,-bac_order,-bac_phylum) %>% 
  dplyr::filter(label %in% groups_to_include) %>%
  group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))



# 4. get together metaxa and idtaxa
rdna<-rbind(idtaxa_occ_target,metaxa_occ_target) %>% select(-n) %>% distinct() %>%
  pivot_wider(values_from = occurrence, values_fill = 0, names_from = label) %>%
  filter(metagenome %in% metagenomes_with_myco)
rdna<- data_frame("type"=c(rep("assembly",length(metagenomes_with_myco)),
                           rep("reads",length(metagenomes_with_myco))), 
                  "metagenome"=c(metagenomes_with_myco,metagenomes_with_myco)) %>% left_join(rdna) %>%
  mutate_if(is.numeric, ~replace_na(.,0))


write.table(rdna,"analysis/03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv",sep="\t",quote = F, row.names = F)


## 5. make a heat map with three levels of detection by mycobiont group. this time use the 330 dataset!

#### add labels: on the level of classes for non-lecanoromycetes, orders for non-lecanorales lecanoromycetes
#### lecanorales split into parmeliaceaea, cladoniaceae, and rest
mtg_info <- read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
mtg_info$label<-mtg_info$order
mtg_info2 <- mtg_info %>% 
  mutate(label_myco=ifelse(class != "Lecanoromycetes", class,
                           ifelse(order != "Lecanorales",order,
                                  ifelse(family %in% c("Parmeliaceae","Cladoniaceae"), family,"Lecanorales_rest")))) %>%
  select(Run,label_myco,photobiont)

####combine mag and rdna tables + add labels
mag_presence$type<-"MAG"
mag_presence$Ulvophyceae<-0
presence<-rbind(mag_presence,rdna) %>%
  pivot_longer(-c(metagenome,type),names_to="group",values_to = "occ") %>%
  left_join(mtg_info2,by=c("metagenome"="Run")) 

##change one label
presence$label_myco[presence$label_myco=="Lecanoromycetes ins. ced."]<-"Lecanorales_rest"

##filter out

####summarize separately by mycobiont and photobiont, calculate the occurrence % and combine together
myco <- presence %>% group_by(type,group,label_myco) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=label_myco,grouping="myco") %>% select(-label_myco)
photo <- presence %>% group_by(type,group,photobiont) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=photobiont,grouping="photo") %>% select(-photobiont)
#calculate n, not occ rate
photo2 <- presence %>% group_by(type,group,photobiont) %>% summarize(found_in=sum(occ),total_mtg=n()) %>%
  mutate(label=photobiont,grouping="photo") %>% select(-photobiont)

summarized<-rbind(myco,photo) %>% mutate(symbiont_status=ifelse(group %in% c("Trebouxiophyceae","Ulvophyceae","Cyanobacteria"),"established","non-established"))


####rename and reorder factors
summarized$type<-factor(summarized$type,levels=c("MAG","assembly","reads"))
summarized$group<-factor(summarized$group,
                       levels=c("Trebouxiophyceae","Ulvophyceae","Cyanobacteria","Acetobacteraceae","Acidobacteriaceae","Beijerinckiaceae","Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes"))

#summarized$group<-factor(summarized$group,
 #                        levels=c("Trebouxiaceae","Coccomyxa-Elliptochloris clade","Ulvophyceae","Cyanobacteria","Acetobacteraceae","Acidobacteriaceae","Beijerinckiaceae","Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes"))

summarized$grouping<-factor(summarized$grouping,levels=c("myco","photo"))

summarized$label[summarized$label=="cyano"]<-"Cyanobacterial"
summarized$label[summarized$label=="trebouxioid"]<-"Trebouxioid"
summarized$label[summarized$label=="trentepohlioid"]<-"Trentepohlioid (Ulvophyceae)"
summarized$label[summarized$label=="trebouxioid_cyano"]<-"Trebouxioid + Cyanobacterial"
summarized$label<-factor(summarized$label,
                            levels=c("Trebouxioid + Cyanobacterial","Cyanobacterial","Trentepohlioid (Ulvophyceae)","Trebouxioid",
                                     "Lecanorales_rest","Cladoniaceae","Parmeliaceae",
                                     "Teloschistales","Caliciales","Leprocaulales",
                                     "Peltigerales","Lecideales","Rhizocarpales",
                                     "Baeomycetales","Gyalectales","Ostropales",
                                     "Pertusariales","Schaereriales" ,"Ostropomycetidae ins. ced.",
                                     "Sarrameanales", "Umbilicariales","Acarosporales",
                                     "Eurotiomycetes","Dothideomycetes","Arthoniomycetes",
                                     "Lichinomycetes"))



####visualize
###filter out groups with <4 samples representd in the 330 dataset

mags<-ggplot(summarized %>% filter(type=="MAG",n>3,!(is.na(label))), 
       aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("MAGs")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~symbiont_status,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),strip.text = element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

assemblies<-ggplot(summarized %>% filter(type=="assembly",n>3,!(is.na(label))), 
             aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Assemblies")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~symbiont_status,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),strip.text = element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

reads<-ggplot(summarized %>% filter(type=="reads",n>3,!(is.na(label))), 
                   aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Reads")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~symbiont_status,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),strip.text = element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

mags+ plot_spacer() +assemblies+ plot_spacer() +reads+plot_layout(widths = c(1,-0.15,1,-0.15,1))

ggsave("results/multilevel_heatmap.png",width=9,height=5,bg="white")
ggsave("results/multilevel_heatmap.pdf",width=9,height=5,bg="white")
ggsave("results/multilevel_heatmap.svg",width=9,height=5,bg="white")


## 6. similar heatmap for ind bacterial genera

####select top bacterial genera, which are:
##among the 13 genera selected for the funcational annotation
##are present in the idtaxa database (which didn't have some of the genera like CAHJXG01)
selected_genera<-mags_role %>% filter(!is.na(bac_genus2)) %>% group_by(bac_genus2) %>%
  summarize(n=n()) %>% arrange(desc(n)) %>% filter(bac_genus2!="unclassified Acetobacteraceae")
selected_genera<-selected_genera$bac_genus2[1:12]
selected_genera<-selected_genera[selected_genera %in% idtaxa$bac_genus2]

####get their occurrence in mags
mag_presence_genera <- mags_role %>% dplyr::select(metagenome,breadth,bac_genus2) %>%
  dplyr::filter(bac_genus2 %in% selected_genera) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,bac_genus2 ) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = bac_genus2 )
####add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence_genera <- data_frame("metagenome"=metagenomes_with_myco) %>% left_join(mag_presence_genera) %>%
  mutate_if(is.numeric, ~replace_na(.,0))
mag_presence_genera$type<-"MAG"

####in rdna
####get genus_level assignment
bac_genus<-plyr::ldply(l_tmp, rbind)[8]
idtaxa$bac_genus2<-bac_genus[,1]

metagenomes_with_mags<-mags_role$metagenome %>% unique()

idtaxa_occ_genera<-idtaxa %>% dplyr::select(-assignment,-file,-bac_order,-bac_phylum,-label) %>% 
  dplyr::filter(bac_genus2 %in% selected_genera) %>%
  group_by(type,metagenome,bac_genus2) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0)) %>%
  pivot_wider(-n,values_from = occurrence, values_fill = 0, names_from = bac_genus2) %>%
  filter(metagenome %in% metagenomes_with_mags)
idtaxa_occ_genera<- data_frame("type"=c(rep("assembly",length(metagenomes_with_myco)),
                                        rep("reads",length(metagenomes_with_myco))), 
                               "metagenome"=c(metagenomes_with_myco,metagenomes_with_myco)) %>% left_join(idtaxa_occ_genera) %>%
  mutate_if(is.numeric, ~replace_na(.,0))

####combine
presence_genera<-rbind(mag_presence_genera,idtaxa_occ_genera) %>%
  pivot_longer(-c(metagenome,type),names_to="group",values_to = "occ") %>%
  left_join(mtg_info2,by=c("metagenome"="Run"))

##change one label
presence_genera$label_myco[presence_genera$label_myco=="Lecanoromycetes ins. ced."]<-"Lecanorales_rest"

####summarize separately by mycobiont and photobiont, calculate the occurrence % and combine together
myco_genera <- presence_genera %>% group_by(type,group,label_myco) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=label_myco,grouping="myco") %>% select(-label_myco)
photo_genera <- presence_genera %>% group_by(type,group,photobiont) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=photobiont,grouping="photo") %>% select(-photobiont)

summarized_genera<-rbind(myco_genera,photo_genera) %>% filter(n>3)

####rename and reorder factors
summarized_genera$type<-factor(summarized_genera$type,levels=c("MAG","assembly","reads"))

summarized_genera$grouping<-factor(summarized_genera$grouping,levels=c("myco","photo"))

summarized_genera$label[summarized_genera$label=="cyano"]<-"Cyanobacterial"
summarized_genera$label[summarized_genera$label=="trebouxioid"]<-"Trebouxioid"
summarized_genera$label[summarized_genera$label=="trentepohlioid"]<-"Trentepohlioid (Ulvophyceae)"
summarized_genera$label[summarized_genera$label=="trebouxioid_cyano"]<-"Trebouxioid + Cyanobacterial"
summarized_genera$label<-factor(summarized_genera$label,
                         levels=c("Trebouxioid + Cyanobacterial","Cyanobacterial","Trentepohlioid (Ulvophyceae)","Trebouxioid",
                                  "Lecanorales_rest","Cladoniaceae","Parmeliaceae",
                                  "Teloschistales","Caliciales","Leprocaulales",
                                  "Peltigerales","Lecideales","Rhizocarpales",
                                  "Baeomycetales","Gyalectales","Ostropales",
                                  "Pertusariales","Schaereriales" ,"Ostropomycetidae ins. ced.",
                                  "Sarrameanales", "Umbilicariales","Acarosporales",
                                  "Eurotiomycetes","Dothideomycetes","Arthoniomycetes",
                                  "Lichinomycetes"))



####visualize
###filter out groups with <4 samples


mags_genera<-ggplot(summarized_genera %>% filter(type=="MAG",n>=4,!(is.na(label))), 
             aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("MAGs")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

assemblies_genera<-ggplot(summarized_genera %>% filter(type=="assembly",n>=4,!(is.na(label))), 
                   aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Assemblies")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

reads_genera<-ggplot(summarized_genera %>% filter(type=="reads",n>=4,!(is.na(label))), 
              aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Reads")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  scale_fill_gradient( low = "white", high = "#3A3A98")+
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

mags_genera+ plot_spacer() +assemblies_genera+ plot_spacer() +reads_genera+plot_layout(widths = c(1,-0.15,1,-0.15,1))

ggsave("results/multilevel_heatmap_genera.png",width=9,height=5,bg="white")
ggsave("results/multilevel_heatmap_genera.svg",width=9,height=5,bg="white")
ggsave("results/multilevel_heatmap_genera.pdf",width=9,height=5,bg="white")

## 7. Visualize less frequent genera

####select top bacterial genera, which are:
##among the 13 genera selected for the funcational annotation
##are present in the idtaxa database (which didn't have some of the genera like CAHJXG01)
selected_genera<-c("Sphingomonas_G","Sphingomonas_I","Sphingomonas_N","Methylobacterium",
                    "Enterovirga","HMF7410","Hymenobacter")


####get their occurrence in mags
mag_presence_genera <- mags_role %>% dplyr::select(metagenome,breadth,bac_genus2) %>%
  dplyr::filter(bac_genus2 %in% selected_genera) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,bac_genus2 ) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = bac_genus2 )
####add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence_genera <- data_frame("metagenome"=metagenomes_with_myco) %>% left_join(mag_presence_genera) %>%
  mutate_if(is.numeric, ~replace_na(.,0))
mag_presence_genera$type<-"MAG"

####in rdna
####get genus_level assignment
bac_genus<-plyr::ldply(l_tmp, rbind)[8]
idtaxa$bac_genus2<-bac_genus[,1]

idtaxa_occ_genera<-idtaxa %>% dplyr::select(-assignment,-file,-bac_order,-bac_phylum,-label) %>% 
  dplyr::filter(bac_genus2 %in% selected_genera) %>%
  group_by(type,metagenome,bac_genus2) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0)) %>%
  pivot_wider(-n,values_from = occurrence, values_fill = 0, names_from = bac_genus2) %>%
  filter(metagenome %in% metagenomes_with_mags)
idtaxa_occ_genera<- data_frame("type"=c(rep("assembly",length(metagenomes_with_myco)),
                                        rep("reads",length(metagenomes_with_myco))), 
                               "metagenome"=c(metagenomes_with_myco,metagenomes_with_myco)) %>% left_join(idtaxa_occ_genera) %>%
  mutate_if(is.numeric, ~replace_na(.,0))

####combine
presence_genera<-rbind(mag_presence_genera,idtaxa_occ_genera) %>%
  pivot_longer(-c(metagenome,type),names_to="group",values_to = "occ") %>%
  left_join(mtg_info2,by=c("metagenome"="Run"))

##change one label
presence_genera$label_myco[presence_genera$label_myco=="Lecanoromycetes ins. ced."]<-"Lecanorales_rest"

####summarize separately by mycobiont and photobiont, calculate the occurrence % and combine together
myco_genera <- presence_genera %>% group_by(type,group,label_myco) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=label_myco,grouping="myco") %>% select(-label_myco)
photo_genera <- presence_genera %>% group_by(type,group,photobiont) %>% summarize(occ_rate=round(mean(occ),2),n=n()) %>%
  mutate(label=photobiont,grouping="photo") %>% select(-photobiont)

summarized_genera<-rbind(myco_genera,photo_genera) %>% filter(n>3)

####rename and reorder factors
summarized_genera$type<-factor(summarized_genera$type,levels=c("MAG","assembly","reads"))

summarized_genera$grouping<-factor(summarized_genera$grouping,levels=c("myco","photo"))

summarized_genera$label[summarized_genera$label=="cyano"]<-"Cyanobacterial"
summarized_genera$label[summarized_genera$label=="trebouxioid"]<-"Trebouxioid"
summarized_genera$label[summarized_genera$label=="trentepohlioid"]<-"Trentepohlioid (Ulvophyceae)"
summarized_genera$label[summarized_genera$label=="trebouxioid_cyano"]<-"Trebouxioid + Cyanobacterial"
summarized_genera$label<-factor(summarized_genera$label,
                                levels=c("Trebouxioid + Cyanobacterial","Cyanobacterial","Trentepohlioid (Ulvophyceae)","Trebouxioid",
                                         "Lecanorales_rest","Cladoniaceae","Parmeliaceae",
                                         "Teloschistales","Caliciales","Leprocaulales",
                                         "Peltigerales","Lecideales","Rhizocarpales",
                                         "Baeomycetales","Gyalectales","Ostropales",
                                         "Pertusariales","Schaereriales" ,"Ostropomycetidae ins. ced.",
                                         "Sarrameanales", "Umbilicariales","Acarosporales",
                                         "Eurotiomycetes","Dothideomycetes","Arthoniomycetes",
                                         "Lichinomycetes"))



####visualize
###filter out groups with <5 samples


mags_genera<-ggplot(summarized_genera %>% filter(type=="MAG",n>3,!(is.na(label))), 
                    aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("MAGs")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

assemblies_genera<-ggplot(summarized_genera %>% filter(type=="assembly",n>3,!(is.na(label))), 
                          aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Assemblies")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),
        axis.text.y=element_blank(),legend.position="none",
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

reads_genera<-ggplot(summarized_genera %>% filter(type=="reads",n>3,!(is.na(label))), 
                     aes(x = group, y = label, fill = occ_rate)) +
  geom_tile(color = "black") + xlab("")+ylab("")+ggtitle("Reads")+
  geom_text(aes(label = occ_rate), color = "white", size = 2) +
  facet_grid(grouping~.,space = "free",scales="free")+theme_minimal()+
  theme(axis.text.x = element_text(angle=90),strip.text = element_blank(),
        axis.text.y=element_blank(),
        panel.spacing.x = unit(0, "cm", data = NULL), panel.spacing.y=unit(0.8, "cm", data = NULL),
        plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0,0,0,0), "cm"))

mags_genera+ plot_spacer() +assemblies_genera+ plot_spacer() +reads_genera+plot_layout(widths = c(1,-0.15,1,-0.15,1))
mags_genera+ assemblies_genera +reads_genera



## 8. Misc: calculate the difference in detection rate based on diff. methodss
###for high-level groups
diff <- summarized %>% pivot_wider(names_from = type,values_from = occ_rate) %>%
  mutate(assembly_mag = assembly - MAG, reads_assembly = reads - assembly,
         reads_mag = reads - MAG)

#how many are in negative? only three, all in very small groups
diff %>% filter(assembly_mag<0) %>% select(group,n,label,MAG,assembly,assembly_mag)
diff %>% filter(reads_assembly<0) %>% select(group,n,label,assembly,reads, reads_assembly)
diff %>% filter(reads_mag<0) %>% select(group,n,label,MAG,reads, reads_mag)

###for genera
diff2 <- summarized_genera %>% pivot_wider(names_from = type,values_from = occ_rate) %>%
  mutate(assembly_mag = assembly - MAG, reads_assembly = reads - assembly,
         reads_mag = reads - MAG)

#how many are in negative? only three, all in very small groups
diff2 %>% filter(assembly_mag<0) %>% select(group,n,label,MAG,assembly,assembly_mag)
diff2 %>% filter(reads_assembly<0) %>% select(group,n,label,assembly,reads, reads_assembly)
diff2 %>% filter(reads_mag<0) %>% select(group,n,label,MAG,reads, reads_mag)



### try with level 7 predictions and split green algae into groups
# 1. mag presence with algal MAGs split by genus
mag_presence <- mags_role %>% dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% groups_to_include) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label)

##add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence <- data_frame("metagenome"=metagenomes_with_myco) %>% left_join(mag_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0))



# 2. get info from metaxa for 18S of CYpho, Tremella, Ulvophyceae, and Trebouxia. screening of assemblies and reads
matrix_level7<-read.delim("analysis/03_metagenome_reanalysis/metaxa_level_7_combined.txt")
matrix_level7_long<-matrix_level7 %>% pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% mutate(occurrence=ifelse(abundance>0,1,0))

##summarize by group, reduce the table only to target groups 
metaxa_occ_target<-matrix_level7_long %>% dplyr::filter(occurrence>0, grepl( "Cystobasidiomycetes", Taxa) |
                                                          grepl( "Tremellomycetes", Taxa) |
                                                          grepl( "Trebouxia", Taxa) |  grepl( "Asterochloris", Taxa) | grepl( "Coccomyxa", Taxa) | grepl( "Elliptochloris", Taxa) |  grepl("Dictyochloropsis", Taxa) | 
                                                          grepl( "Ulvophyc", Taxa)) %>%
  mutate(label=ifelse(grepl("Cystobasidiomycetes",Taxa),"Cystobasidiomycetes",
                      ifelse(grepl("Beijerinckiaceae",Taxa), "Beijerinckiaceae",
                             ifelse(grepl("Tremellomycetes", Taxa),"Tremellomycetes",
                                    ifelse(grepl("Trebouxia", Taxa)|grepl( "Asterochloris", Taxa)|grepl("Dictyochloropsis", Taxa),"Trebouxiaceae",
                                           ifelse(grepl( "Ulvophyc", Taxa),"Ulvophyceae",
                                                  ifelse(grepl( "Coccomyxa", Taxa) | grepl( "Elliptochloris", Taxa),"Coccomyxa-Elliptochloris clade",NA))))))) %>%
  dplyr::select(-Taxa,-abundance) %>% group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))

### add CHloroplast predictions (couldn't be summarized by metaxa itself because of ranking problems, did it maunaly with grep)
chlorop<-read.delim("analysis/03_metagenome_reanalysis/chloroplast_taxonomy.txt",header=F)
chlorop_assembly<-chlorop %>% filter(grepl("assembly",V1)) %>%
  separate(V1,into=c("type","rest"),sep="_") %>% mutate(metagenome=gsub("\\..*","",rest)) %>%
  filter( grepl( "Trebouxia", V2) |  grepl( "Asterochloris", V2) | grepl( "Coccomyxa", V2) | grepl( "Elliptochloris", V2) |  grepl("Dictyochloropsis", V2)) %>%
  mutate( label=ifelse(grepl("Trebouxia", V2)|grepl( "Asterochloris", V2)|grepl("Dictyochloropsis", V2),"Trebouxiaceae",
                                           ifelse(grepl( "Ulvophyc", V2),"Ulvophyceae",
                                                  ifelse(grepl( "Coccomyxa", V2) | grepl( "Elliptochloris", V2),"Coccomyxa-Elliptochloris clade",NA)))) %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>% mutate(occurrence=1)
chlorop_assembly$type<-"assembly"

chlorop_reads<-chlorop %>% filter(!grepl("assembly",V1)) %>%
  mutate(metagenome=gsub("\\..*","",V1)) %>%
  filter( grepl( "Trebouxia", V2) |  grepl( "Asterochloris", V2) | grepl( "Coccomyxa", V2) | grepl( "Elliptochloris", V2) |  grepl("Dictyochloropsis", V2)) %>%
  mutate( label=ifelse(grepl("Trebouxia", V2)|grepl( "Asterochloris", V2)|grepl("Dictyochloropsis", V2),"Trebouxiaceae",
                       ifelse(grepl( "Ulvophyc", V2),"Ulvophyceae",
                              ifelse(grepl( "Coccomyxa", V2) | grepl( "Elliptochloris", V2),"Coccomyxa-Elliptochloris clade",NA)))) %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>% mutate(occurrence=1)
chlorop_reads$type<-"reads"

  
metaxa_occ_target<-rbind(metaxa_occ_target,chlorop_reads,chlorop_assembly)


##try removing chloroplast genes
metaxa_occ_target<-matrix_level7_long %>% dplyr::filter(occurrence>0, grepl( "Cystobasidiomycetes", Taxa) |
                                                          grepl( "Tremellomycetes", Taxa) |
                                                          grepl( "Trebouxia", Taxa) |  grepl( "Asterochloris", Taxa) | grepl( "Coccomyxa", Taxa) | grepl( "Elliptochloris", Taxa) |  grepl("Dictyochloropsis", Taxa) | 
                                                          grepl( "Ulvophyc", Taxa)) %>%
  mutate(label=ifelse(grepl("Cystobasidiomycetes",Taxa),"Cystobasidiomycetes",
                      ifelse(grepl("Beijerinckiaceae",Taxa), "Beijerinckiaceae",
                             ifelse(grepl("Tremellomycetes", Taxa),"Tremellomycetes",
                                    ifelse(grepl("Trebouxia", Taxa)|grepl( "Asterochloris", Taxa)|grepl("Dictyochloropsis", Taxa),"Trebouxiaceae",
                                           ifelse(grepl( "Ulvophyc", Taxa),"Ulvophyceae",
                                                  ifelse(grepl( "Coccomyxa", Taxa) | grepl( "Elliptochloris", Taxa),"Coccomyxa-Elliptochloris clade",NA))))))) %>%
  dplyr::select(-Taxa,-abundance) %>% group_by(type,metagenome,label) %>%
  summarize(n=n()) %>% mutate(occurrence=ifelse(n>0,1,0))
