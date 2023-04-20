#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(patchwork)
library(scales)
source('code/utils.R')

#load in data
##read table with metadata
metatable<-read.delim2('analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt')

##read table with bais pair counts
bp<-read.delim('analysis/03_metagenome_reanalysis/bp_report.txt',header=F)
colnames(bp)<-c("Run","bp")

##read results of search in raw reads
reads<-read.delim('analysis/03_metagenome_reanalysis/metaxa_report_reads.txt',header=F)
colnames(reads)<-c("Run","yeast","hit_reads")

##read results of search in assembly
assemblies<-read.delim('analysis/03_metagenome_reanalysis/metaxa_report_assembly.txt',header=F)
colnames(assemblies)<-c("Run","yeast","hit_assembly")


#combine data
##join the tables
df<-left_join(bp,reads) %>% left_join(assemblies) %>% 
  mutate(hit_reads_abs=ifelse(hit_reads=="Y","present","absent"),hit_assembly_abs=ifelse(hit_assembly=="Y","present","absent")) %>%
  mutate(presence=ifelse(hit_assembly_abs=="present","present",ifelse(hit_reads_abs=="present","present","absent"))) %>%
  left_join(metatable) %>% select(Run, Source,bp, yeast,presence,Lichen.metagenomes,architecture,class) %>%
  pivot_wider(names_from = yeast, values_from = presence)

##dot histogram
dothist<-ggplot(df, aes(bp,fill=class)) + geom_dotplot(binwidth=250000000,stackgroups = TRUE, binpositions="all", method='histodot') +
  scale_y_continuous(name = "Number of Metagenomes", breaks = NULL)+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("<250 Mbp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
   theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_dothist_by_source.png",device="png",height = 5,width=10)

