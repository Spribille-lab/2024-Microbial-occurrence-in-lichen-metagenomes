# Geographic maps showing detection of symbionts by mag, assembly, reads
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)

# 1. Load geographic info
coord<-read.delim("analysis/03_metagenome_reanalysis/locations_final.txt")


# 2. Get occurrences by MAGs for the selected groups
## load the big table
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_family2=ifelse(bac_family=="Unknown",paste(bac_order," fam.",sep=""),bac_family))
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family2," gen. sp.",sep=""),bac_genus))
mags_role$bac_genus2[mags_role$bac_genus2=="Acetobacteraceae gen. sp." & !is.na(mags_role$bac_genus2)]<-"unclassified Acetobacteraceae"

##add labels
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Trebouxiophyceae"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cystobasidiomycetes"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremellomycetes"

## prepare matrix for the selected groups: cystobasidiomycetes and tremellomycetes and bactera on the level of families
groups_to_include<-c("Trebouxiophyceae","Acidobacteriaceae",
                     "Beijerinckiaceae","Acetobacteraceae",
                     "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes")

mag_presence_fam <- mags_role %>% dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% groups_to_include) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() 

## do the same for selected bacterial genera
selected_genera <- c("Lichenihabitans","CAHJXG01","unclassified Acetobacteraceae","EB88",
                     "CAIMSN01","CAHJWO01","Terriglobus","Nostoc")
mag_presence_genera <- mags_role %>% mutate(label=bac_genus2) %>%
  dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% selected_genera) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() 

## combine
mag_presence<-rbind(mag_presence_fam,mag_presence_genera)
mag_presence$type="genome"

# 3. Get occurrences by rRNA for the selected groups
## already have the matrix for the bac_family and eukaryotes
rdna_fam<-read.delim("analysis/03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv")
rdna_fam <- rdna_fam %>% pivot_longer(-c(type,metagenome), names_to = "label",values_to = "n") %>% filter(n>0) %>% select(-n)


## prepare the matrix for selected genera
idtaxa_assemblies<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt")
idtaxa_assemblies$type<-"assembly"
idtaxa_reads<-read.delim("analysis/03_metagenome_reanalysis/idtaxa_reads.txt")
idtaxa_reads$type<-"reads"
idtaxa<-rbind(idtaxa_assemblies,idtaxa_reads) %>% filter(metagenome %in% metagenomes_with_mags) #added filtering step to make sure that only intended metagenomes are counted, i.e. the metagenomes with a mycobiont mag and those that were not excluded as misidentified

##get genus_level assignment
l_tmp<-str_split(idtaxa$assignment,";") 
bac_genus<-plyr::ldply(l_tmp, rbind)[8]
idtaxa$label<-bac_genus[,1]
rdna_genera<-idtaxa %>% filter(label %in% selected_genera) %>% select(metagenome,type,label) %>% distinct()

## 4. Combine tables
df<-rbind(rdna_genera,rdna_fam,mag_presence)
df$presence<-1
df2<-df %>% pivot_wider(names_from = label,values_from = presence,values_fill = 0) %>%
  pivot_longer(-c(type,metagenome),names_to = "label",values_to="presence") %>%
  pivot_wider(names_from = type,values_from = presence,values_fill = 0) %>%
  pivot_longer(-c(label,metagenome),names_to = "type",values_to="presence")

## add geography
df2<-df2 %>% left_join(coord,by=c("metagenome"="run_accession"))
df2$metagenome<-as.factor(df2$metagenome)
df2$data_type<-factor(df2$type,levels=c("genome","assembly","reads"))
df2$label<-as.factor(df2$label)

# 5. make a map
world <- ne_countries(scale = "medium", returnclass = "sf")

##all lineages on one turned out to be too much for one figure, so split into two
df3<-df2 %>% filter(label %in% groups_to_include | label=="Lichenihabitans")

ggplot(data = world) +
  geom_sf(fill= "antiquewhite",color="lightgrey",size=0.2)+
  coord_sf(expand = FALSE)+
  geom_point(data = df3, aes(x = as.numeric(long), y = as.numeric(lat),alpha=presence), color = "brown", size = 2)+ 
  theme_minimal()+
  scale_alpha(range = c(0, 0.3))+
  facet_grid(label~data_type)+
  theme(axis.title=element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "aliceblue"),
        strip.text = element_text(angle = 0))
  
ggsave("analysis/03_metagenome_reanalysis/map_symbionts.pdf",width=7.8,height=10)
ggsave("analysis/03_metagenome_reanalysis/map_symbionts.png",width=7.8,height=10)
