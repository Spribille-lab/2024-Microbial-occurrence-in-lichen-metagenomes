#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
require(scales)


##define functions
library(scales)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

    
    
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mtg_info<-read.delim("results/tables/all_metagenome_reanalysis.txt")

#add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))
mags_role<-mags_role %>%mutate(bac_family2=ifelse(bac_family=="Unknown",paste(bac_order," fam. gen. sp.",sep=""),bac_genus))

mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremella"

mags_role$label_genus<-mags_role$bac_genus2
mags_role$label_genus[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label_genus[mags_role$confirmed_role=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label_genus[mags_role$confirmed_role=="cypho"]<-"Cyphobasidium"
mags_role$label_genus[mags_role$confirmed_role=="tremella"]<-"Tremella"

mags_role$label_genus_level<-mags_role$bac_family
mags_role$label_genus_level[mags_role$confirmed_role=="photobiont_chloro"]<-"Eukaryotes"
mags_role$label_genus_level[mags_role$confirmed_role=="cypho"]<-"Eukaryotes"
mags_role$label_genus_level[mags_role$confirmed_role=="tremella"]<-"Eukaryotes"
mag_taxonomy<-mags_role %>% dplyr::select(bac_genus2,bac_family2,bac_order,bac_phylum) %>%
  distinct()
  
#filter: remove metagenomes that didn't yeild mycobiont mag
#remove "fungi_other" 
mags_role_filtered<-mags_role %>% filter(!(confirmed_role %in% c("fungi_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Mycobiont"])


#add column for the depth of coverage of the mycobiont and for the coverage relative to the myco_cov
mags_role_filtered<-mags_role_filtered %>% group_by(metagenome) %>% 
  mutate(myco_cov = ifelse(confirmed_role=="mycobiont",depth_cov,NA)) %>%
  fill(myco_cov, .direction = 'updown') %>% mutate(relative_cov = depth_cov/myco_cov) %>%
  ungroup()


#plot by familiy
groups_to_include=c("Beijerinckiaceae","Acetobacteraceae","Sphingomonadaceae","Acidobacteriaceae","Nostocaceae","Cyphobasidium","Tremella","Green Algal Photobiont")
ggplot(mags_role_filtered %>% filter(label %in% groups_to_include),aes(x=label,y=relative_cov))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(trans = squish_trans(2, 14, 10),
                     breaks = seq(0, 14, by = 2))+
  geom_hline(yintercept = 1, colour = "red")+
  xlab("")+ylab("Depth of coverage relative \n to the mycobiont MAG")


#plot by genus
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
annotated_mags<-annotated_mags %>% left_join(mags_role %>% dplyr::select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())

ggplot(mags_role_filtered %>% filter(label_genus %in% annotated_mags$bac_genus | label_genus_level=="Eukaryotes"),aes(x=label_genus,y=relative_cov,fill=label_genus_level,col=label_genus_level))+
  geom_boxplot()+
  scale_y_continuous(trans = squish_trans(1.5, 14, 40),
                     breaks = c(0.25,0.5,0.75,1,1.25,seq(2, 14, by = 4)))+
  geom_hline(yintercept = 1, colour = "red")+
  facet_grid(.~label_genus_level,space="free",scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=18),
        strip.text.x = element_blank(),legend.title = element_blank())+
  xlab("")+ylab("Depth of coverage relative to \n the main fungal symbiont MAG")

ggsave("results/figures/relative_cov_boxplot_genus.pdf",bg="white",width=18,height=9)

#how many mags had relative cov >1?
mags_role_filtered %>% filter(relative_cov>1) %>%
  group_by(bac_order) %>% summarize(n=n())

mags_role_filtered %>% filter(relative_cov>1) %>%
  group_by(label) %>% summarize(n=n())


high_cov<-mags_role_filtered %>% filter(label!="Green Algal Photobiont" & label!="Mycobiont" & relative_cov>0.5)

1-nrow(high_cov)/nrow(mags_role_filtered %>% filter(label!="Green Algal Photobiont" & label!="Mycobiont"))

#is nostoc mostly >1 or <1
mags_role_filtered %>% filter(bac_genus=="Nostoc") %>% 
  mutate(binary=ifelse(relative_cov>1,"more","les")) %>%
  group_by(binary) %>% summarize(n=n())


##coverage relative to summed bacterial coverage
mags_role_filtered2<-mags_role_filtered %>% group_by(metagenome) %>% 
  filter(lineage_broad %in% c("bacteria_other","Cyanobacteria")) %>%
  mutate(total_bac_cov = sum(depth_cov)) %>%
  fill(total_bac_cov, .direction = 'updown') %>% mutate(relative_cov_bac = depth_cov/total_bac_cov) %>%
  ungroup()
ggplot(mags_role_filtered2 %>% filter(bac_order %in% orders_to_include),aes(x=bac_genus,y=relative_cov_bac))+
  geom_violin()+geom_jitter(width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~bac_order,scales = "free_x")+xlab("")+ylab("Depth of coverage relative \n to the mycobiont MAG")
ggsave("results/figures/bac_cov_relative_to_all_bac.png",bg="white",width=18,height=7)

## average relative coverage by genus
mags_role_filtered %>%  filter(label_genus %in% annotated_mags$bac_genus | label_genus_level=="Eukaryotes") %>%
  group_by(label_genus) %>% summarize(avg_rel_cov=mean(relative_cov),median=median(relative_cov),max=max(relative_cov),min=min(relative_cov))

##distribution of relative coverage of algal mags
ggplot(mags_role_filtered %>% filter(label_genus=="Green Algal Photobiont"))+
  geom_histogram(aes(x=relative_cov),binwidth = 0.05)+
  scale_x_continuous(breaks=seq(0,1.5,by=0.1))

ggplot(mags_role_filtered %>% filter(label_genus=="Nostoc"))+
  geom_histogram(aes(x=relative_cov))+
  scale_x_continuous(breaks=seq(0,20,by=1))

ggplot(mags_role_filtered %>% filter(label_genus=="CAHJXG01"))+
  geom_histogram(aes(x=relative_cov),binwidth = 0.05)+
  scale_x_continuous(breaks=seq(0,1.5,by=0.1))

ggplot(mags_role_filtered %>% filter(label_genus=="Lichenihabitans"))+
  geom_histogram(aes(x=relative_cov),binwidth = 0.05)+
  scale_x_continuous(breaks=seq(0,5,by=0.2))


###make table for export
rel_cov_df<-mags_role_filtered %>% group_by(label_genus) %>% filter(label_genus!="Mycobiont") %>%
  summarize(median_relative_cov=median(relative_cov)) %>% left_join(mag_taxonomy,by=c("label_genus"="bac_genus2"))
write.table(rel_cov_df,"results/tables/median_relative_coverage_by_genus.tsv",sep="\t",quote = F, row.names = F)

