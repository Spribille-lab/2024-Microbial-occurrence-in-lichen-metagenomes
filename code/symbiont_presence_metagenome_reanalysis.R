#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(patchwork)
library(scales)
library(R.utils)
source('code/utils.R')

#load in data
##read table with metadata
metatable<-read.delim2('analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt')

##read table with bais pair counts
bp<-read.delim('analysis/03_metagenome_reanalysis/bp_report.txt',header=F)
colnames(bp)<-c("Run","bp")

##read results of search in raw reads
reads<-read.delim('analysis/03_metagenome_reanalysis/metaxa_report_reads.txt',header=F)
colnames(reads)<-c("Run","symbiont","hit_reads")

##read results of search in assembly
assemblies<-read.delim('analysis/03_metagenome_reanalysis/metaxa_report_assembly.txt',header=F)
colnames(assemblies)<-c("Run","symbiont","hit_assembly")


#combine data
##join the tables
df<-left_join(bp,reads) %>% left_join(assemblies) %>% 
  mutate(hit_reads_abs=ifelse(hit_reads=="Y","present","absent"),hit_assembly_abs=ifelse(hit_assembly=="Y","present","absent")) %>%
  mutate(presence=ifelse(hit_assembly_abs=="present","present in assembly and reads",ifelse(hit_reads_abs=="present","present only in reads","absent"))) %>%
  left_join(metatable) %>% select(Run, bp, symbiont,presence,Lichen.metagenomes,architecture)
  
df2<-df %>% pivot_wider(names_from = symbiont, values_from = presence)
df2$seq_depth<- hsize(df2$bp,standard="SI")

df2<-df2 %>% select(-bp)

#save the table
write.csv(df2,"results/tables/symbiont_presence_summary.csv",row.names = F,quote = F)

## vis
df$symbiont<-factor(df$symbiont, levels = c("cypho", "tremella","trebouxia"), 
                  labels = c("Cyphobasidium", "Tremella","Trebouxia"))
df$presence<-factor(df$presence, levels = c("present in assembly and reads", "present only in reads","absent"))
dothist<-ggplot(df, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(symbiont~., ncol=1)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist.png",device="png",height = 8,width=10)


