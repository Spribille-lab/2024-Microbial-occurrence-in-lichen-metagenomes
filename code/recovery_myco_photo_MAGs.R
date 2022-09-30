#setwd("~/Documents/coverage")
library(tidyverse)
library(stringr)
library(patchwork)

# identify mycobiont and photobiont recovery as a function of depth
# 1. based on the recovery of genomes from the read alignment data (i.e. aligning reads from all metagenomes to all MAGs). A MAG counts as present, if it has breadth of coverage >50%, according to BWA
### Saved the figure as `results/figures/myco_photobiont_vs_depth.png`

##load data
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")

##select only occurrences of myco- and photobionts
mag_recovery<-mags_role %>% filter(confirmed_role %in% c("mycobiont","mycobiont_missassigned","photobiont_chloro","photobiont_cyano")) %>%
  select(metagenome,confirmed_role) %>% 
  mutate(role_type=ifelse(confirmed_role %in% c("mycobiont","mycobiont_missassigned"),"mycobiont","photobiont")) %>%
  select(-confirmed_role) %>% group_by(metagenome,role_type) %>%
  summarize(n=n()) %>% pivot_wider(names_from = role_type,values_from = n, values_fill=0)


##add sequencing depth information 
depth<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F)
colnames(depth)<-c("metagenome","depth")
mag_recovery<-depth %>% left_join(mag_recovery)

##in the metagenomes that yilded no mags, replace NAs with 0s
mag_recovery$mycobiont[is.na(mag_recovery$mycobiont)]<-0
mag_recovery$photobiont[is.na(mag_recovery$photobiont)]<-0

##add binary variables
mag_recovery<-mag_recovery %>% mutate(mycobiont_binary=ifelse(mycobiont>0,"detected","not detected"),
                        photobiont_binary=ifelse(photobiont>0,"detected","not detected"),
                        photo_and_myco=ifelse(mycobiont>0&photobiont>0,"detected","not detected"))

##visualizing  recovery
### whether myco/photobiont mag was detected vs coverage: jitter
myco_detection<-ggplot(mag_recovery,aes(y=depth,x=mycobiont_binary))+geom_jitter(width = 0.1,alpha=0.5,col="dark red")+
  theme_minimal()+xlab("Lichen fungal symbiont genome")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))

photo_detection<-ggplot(mag_recovery,aes(y=depth,x=photobiont_binary))+geom_jitter(width = 0.1,alpha=0.5,col="dark green")+
  theme_minimal()+xlab("Photosynthetic symbionts genome")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))

both_detection<-ggplot(mag_recovery,aes(y=depth,x=photo_and_myco))+geom_jitter(width = 0.1,alpha=0.5,col="orange4")+
  theme_minimal()+xlab("LFS & Photosynthetic\nsymbionts genomes")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))

myco_detection+photo_detection+both_detection
ggsave("results/figures/myco_photobiont_vs_depth.png",bg="white",width=10,height=5)

#get numbers
## min depth to get mycbiont mag
min(mag_recovery[mag_recovery$mycobiont>0,]$depth)
#248125314 = 248 Mbp

## min depth to get photobiont mag
min(mag_recovery[mag_recovery$photobiont>0,]$depth)
#168207054 = 168 Mbp

## min depth to get mycbiont AND photobiont mag
min(mag_recovery[mag_recovery$photo_and_myco=="detected",]$depth)

min(mags_in_mtg[mags_in_mtg$photo_and_myco==T,]$depth)
#1157540028 = 1.2 Gbp






#2.Based on recovery of MAGs from de-novo binning (including those pre dereplication). If a genome was present according to bwa alignments, but not recovered as MAG (perhaps due to e.g. low coverage), we didn't count it here
##Load data
###load mag information
mags<-read.delim("analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv")

### add euk MAGs that were removed in dereplication
drep<-read.csv("analysis/05_MAGs/tables/drep_euks/data_tables/Cdb.csv")
drep_qc<-read.csv("analysis/05_MAGs/tables/drep_euks/data_tables/genomeInfo.csv")
drep<- drep %>% left_join(drep_qc) %>% filter(completeness>90)

### add functional assignments from the bwa table
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

## Which MAGs were identified as mycobiont and photobiont(s)?
mags_role2<-mags_role %>% filter(confirmed_role %in% c("mycobiont","mycobiont_missassigned","photobiont_chloro","photobiont_cyano")) %>%
        select(Genome,metagenome,confirmed_role) %>% 
  mutate(role_type=ifelse(confirmed_role %in% c("mycobiont","mycobiont_missassigned"),"mycobiont","photobiont")) %>%
  select(-confirmed_role)

## What was the source of this MAG? Was it recovered from a different metagenome?
mags_role2$source_metagenome<-gsub('[a-zA-Z]+_(.*)_.*_.*', '\\1', mags_role2$Genome)
mags_role2<-mags_role2 %>% mutate(is_mtg_same=ifelse(metagenome==source_metagenome,T,F))

## Which metagenomes' myco and photo mags are from metgaenomes different from itself?
mags_role2 %>% group_by(metagenome,role_type,is_mtg_same) %>% summarise(n=n())


## Process dereplication results
drep$metagenome<-gsub('[a-zA-Z]+_(.*)_.*_.*', '\\1', drep$genome)
### link all MAGs to the winning MAG of their cluster
winning_mags<-drep %>% group_by(secondary_cluster) %>% 
  summarize(winning_mag= genome[genome %in% mags$fasta])
drep <- drep %>% left_join(winning_mags)


## This function links a "foreign" MAG present in a metagenome to its version discarded during dereplication
###if there is a discarded MAG, it returns a dataframe that follows the schema of mags_role2
### if there isn't it returns an empty dataframe

linking_mags<-function(t){ #takes a row of mags_role2 table
test<-drep %>% filter(winning_mag==paste(t[1],".fa",sep='') & metagenome == t[2])
if(nrow(test)>0){
  output<-data.frame("Genome"=test$genome,"metagenome"=test$metagenome,"role_type"=t[3],"source_metagenome"=test$metagenome,"is_mtg_same"=T)
}else{
  output<- data.frame("Genome"=character(),"metagenome"=character(),"role_type"=character(),"source_metagenome"=character(),"is_mtg_same"=logical())
}
  return(output)
}

###applying the function to all MAG-metagenome pairs where is_mtg_same==F

l<-apply(mags_role2 %>% filter(is_mtg_same==F),1,linking_mags)
adding_derep<-do.call(rbind,l)


### adding them to the table
mags_derp_combined<-rbind(mags_role2,adding_derep) %>% filter(is_mtg_same==T) %>%
  select(Genome,metagenome,role_type)

#summarize info by metagenome
mags_in_mtg<-mags_derp_combined%>% group_by(metagenome,role_type) %>% summarise(n=n()) %>% 
  pivot_wider(names_from=role_type, values_from=n, values_fill=0)

#add metagenomes that didn't produce any mags
mtg<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")
no_mag_mtg<-data.frame(mtg$Run[!(mtg$Run %in% mags_derp_combined$metagenome)])
colnames(no_mag_mtg)<-"metagenome"
mags_in_mtg<-mags_in_mtg %>% plyr:::rbind.fill(no_mag_mtg)
mags_in_mtg[is.na(mags_in_mtg)] <- 0

#create extra columns for counting various types of mags
mags_in_mtg<-mags_in_mtg %>% mutate(mycobiont_binary=ifelse(mycobiont>0, T, F),photobiont_binary=ifelse(photobiont>0, T, F),photo_and_myco=ifelse(photobiont>0 & mycobiont>0, T, F))
mags_in_mtg<-left_join(mags_in_mtg,mtg,by=c("metagenome"="Run"))         

#add sequencing depth information 
depth<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F)
colnames(depth)<-c("metagenome","depth")
mags_in_mtg<-left_join(mags_in_mtg,depth)
                       

#visualizing  recovery
## whether myco/photobiont mag was detected vs coverage: jitter
myco_detection<-ggplot(mags_in_mtg,aes(y=depth,x=mycobiont_binary))+geom_jitter(width = 0.1,alpha=0.5,col="dark red")+
  theme_minimal()+xlab("Main fungal symbiont MAG detected")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))
  
photo_detection<-ggplot(mags_in_mtg,aes(y=depth,x=photobiont_binary))+geom_jitter(width = 0.1,alpha=0.5,col="dark green")+
  theme_minimal()+xlab("Photosynthetic symbionts MAG detected")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))

both_detection<-ggplot(mags_in_mtg,aes(y=depth,x=photo_and_myco))+geom_jitter(width = 0.1,alpha=0.5,col="orange4")+
  theme_minimal()+xlab("Main fungal & Photosynthetic\nsymbionts MAGs detected")+ylab("Sequencing depth")+
  scale_y_log10(breaks=c(1000000,100000000,10000000000),labels=c("1 Mbp","100 Mbp","10 Gbp"))

myco_detection+photo_detection+both_detection
ggsave("results/figures/myco_photobiont_vs_depth.png",bg="white",width=10,height=5)


#get numbers
## min depth to get mycbiont mag
min(mags_in_mtg[mags_in_mtg$mycobiont_binary==T,]$depth)
#548060842 = 548 Mbp

## min depth to get photobiont mag
min(mags_in_mtg[mags_in_mtg$photobiont_binary==T,]$depth)
# 567138182 = 567 Mbp

## min depth to get mycbiont AND photobiont mag
min(mags_in_mtg[mags_in_mtg$photo_and_myco==T,]$depth)
# 1959059806 = 2 Gbp

