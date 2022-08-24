#visualizeing fegenie results
##are loci identified as pseudogenes located on contig ends?

#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(RColorBrewer)

##mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% dplyr:::select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())


#read fegenie results
fegenie<-read.csv("analysis/07_annotate_MAGs/fegenie/arkadiy/FeGenie-heatmap-data.csv")
colnames(fegenie)[1]<-"fegenie_class"
fegenie$fegenie_class<-row.names(fegenie)

fegenie<-fegenie %>% pivot_longer(-fegenie_class,values_to="n",names_to="fasta")
fegenie$Genome<-fegenie$fasta %>% str_replace(".faa","")
fegenie <- fegenie %>% select(-fasta) #%>% pivot_wider(names_from = fegenie_class,values_from = n)

data <- fegenie %>% left_join(annotated_mags) 


#visualize
data$bac_family_label<-data$bac_family
data[data$bac_family=="Sphingomonadaceae",]$bac_family_label<-"Sphingo\nmonadaceae"
data$fegenie_class<-str_replace_all(data$fegenie_class,"_"," ")
data$fegenie_class<-str_replace_all(data$fegenie_class,"-",":\n")

ggplot(data %>% filter(n>0),aes(color=bac_family_label)) + 
  geom_point(aes(x=Genome,y=fegenie_class,size=n)) +
  ylab("Number of iron metabosm gene clusters")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=18))+
  scale_color_brewer(palette ="Dark2",type="qual")+
  guides(color = "none") +
  facet_grid(.~bac_family_label,scales = "free",space="free",labeller = label_wrap_gen(width=5))
ggsave("results/figures/fegenie_bubbleplot.png",device="png",height = 7,width=16,bg="white")








