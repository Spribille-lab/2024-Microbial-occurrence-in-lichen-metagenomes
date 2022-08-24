#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(patchwork)
library(scales)
source('code/utils.R')
#read data
##read table with metadata
metatable<-read.delim2('analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt')

##read table with bais pair counts
bp<-read.delim('analysis/03_metagenome_reanalysis/bp_report.txt',header=F)
colnames(bp)<-c("Run","bp")

##read results of search in raw reads
reads<-read.delim('analysis/03_metagenome_reanalysis/metaxa_report.txt',header=F)
colnames(reads)<-c("Run","yeast","hit_reads")

##read results of search in assembly
assemblies<-read.delim('analysis/03_metagenome_reanalysis/blast_report.txt',header=F)
colnames(assemblies)<-c("Run","yeast","hit_assembly")

#join the tables
df<-left_join(bp,reads) %>% left_join(assemblies) %>% 
  mutate(hit_reads_abs=ifelse(hit_reads=="Y","present","absent"),hit_assembly_abs=ifelse(hit_assembly=="Y","present","absent")) %>%
  mutate(presence=ifelse(hit_assembly_abs=="present","present in assembly and reads",ifelse(hit_reads_abs=="present","present in reads","absent"))) %>%
  left_join(metatable)

##summarize
df_summary<- df %>% pivot_longer(c(hit_reads_abs, hit_assembly_abs),names_to='method',values_to="outcome") %>% 
  group_by(yeast,method,outcome) %>%   summarise(n = n()) %>% mutate(freq = n / sum(n)) 




##add lines for the two lichen metagenoms we used
G1_trem<-c("G1",10000000000,"tremella","Y","Y","present","present","present in assembly and reads","Bryoria tortuosa","This paper","Spribille Lab, University of Alberta","first","macro","Parmeliaceae","Lecanorales","Lecanoromycetes")
G1_cypho<-c("G1",10000000000,"cypho","Y","Y","present","present","present in assembly and reads","Bryoria tortuosa","This paper","Spribille Lab, University of Alberta","first","macro","Parmeliaceae","Lecanorales","Lecanoromycetes")

my_mtg<-data.frame(rbind(G1_trem,G1_cypho))
names(my_mtg)<-colnames(df)                   
df2<-rbind(df,my_mtg)
df2$bp<-as.numeric(df2$bp)
df2$presence<-factor(df2$presence,levels=c("present in assembly and reads","present in reads","absent"))


#visuazlise dot-histogram
#all, only separated by yeast
df2$yeast<-factor(df2$yeast, levels = c("cypho", "tremella"), 
                             labels = c("Cyphobasidium", "Tremella"))
dothist<-ggplot(df2, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=1)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist.png",device="png",height = 8,width=10)


dothist_low<-ggplot(df2, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,1000000000,2500000000,5000000000),labels=c("0 bp","1 Gbp","2.5 Gbp","5 Gbp"),limits=c(0,5000000000)) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=2)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist_low
ggsave("results/figures/metagenome_reanalysis_dothist_low_depth.png",device="png",height = 8,width=10)





#parmeliaceae, only separated by yeast
df2_parmeliacea<-df2 %>% filter(family=="Parmeliaceae")

dothist<-ggplot(df2_parmeliacea, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=1)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist_parmeliaceae.png",device="png",height = 8,width=10)

dothist_low<-ggplot(df2_parmeliacea, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,1000000000,2500000000,5000000000),labels=c("0 bp","1 Gbp","2.5 Gbp","5 Gbp"),limits=c(0,5000000000)) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=2)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist_low
ggsave("results/figures/metagenome_reanalysis_dothist_parmeliacea_low_depth.png",device="png",height = 8,width=10)


#lecanorales, only separated by yeast
df2_lecanorales<-df2 %>% filter(order=="Lecanorales")
dothist<-ggplot(df2_lecanorales, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=1)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist_lecanorales.png",device="png",height = 8,width=10)

dothist_low<-ggplot(df2_lecanorales, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,1000000000,2500000000,5000000000),labels=c("0 bp","1 Gbp","2.5 Gbp","5 Gbp"),limits=c(0,5000000000)) +
  geom_dotplot(fill=NA,binwidth=200000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~., ncol=2)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist_low
ggsave("results/figures/metagenome_reanalysis_dothist_lecanorales_low_depth.png",device="png",height = 8,width=10)



#grid plot with X axis as arcitecture and Y as the yeast. Only Lecanorales
df2_binary_arc<- df2 %>% filter(architecture %in% c("macro","crust")) %>% filter(order=="Lecanorales")

dothist<-ggplot(df2_binary_arc, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~architecture, ncol=2)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist_lecanorales_arcitechture.png",device="png",height = 8,width=10)


# plot  Y as the yeast. Only Peltigerales
df2_binary_arc<- df2 %>% filter(architecture %in% c("macro","crust")) %>% filter(order=="Peltigerales")

dothist<-ggplot(df2_binary_arc, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(.~yeast, ncol=1)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist_peltigerales.png",device="png",height = 8,width=10)


#grid plot with X axis as arcitecture and Y as the yeast. 
df2_binary_arc<- df2 %>% filter(architecture %in% c("macro","crust")) 

dothist<-ggplot(df2_binary_arc, aes(bp)) + geom_dotplot(fill=col_pres_abs,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot',aes(alpha=presence)) +
  scale_y_continuous(name = "", breaks = NULL)+
  scale_alpha_discrete(range=c(1,0))+
  xlab("Sequencing depth, bp")+
  scale_x_continuous(breaks = c(0,5000000000,10000000000,15000000000,20000000000,25000000000,30000000000,35000000000),labels=c("0 bp","5 Gbp","10 Gbp","15 Gbp","20 Gbp","25 Gbp","30 Gbp","35 Gbp")) +
  geom_dotplot(fill=NA,binwidth=400000000,stackgroups = TRUE, binpositions="all", method='histodot')+
  facet_wrap(yeast~architecture, ncol=2)+
  theme_minimal() + theme(strip.text.x = element_text(size=12,face = "italic"))
dothist
ggsave("results/figures/metagenome_reanalysis_dothist_arcitechture.png",device="png",height = 8,width=10)


##summarize in regards with arcitechture
df_summary2<- df %>% group_by(architecture,yeast,presence) %>%   dplyr::summarise(n = n()) %>% dplyr::mutate(freq = 100*n / sum(n)) %>%
 select(-n) %>%filter(architecture!="Unknown") %>% pivot_wider(names_from = presence,values_from=freq )
write.table(df_summary2,"results/tables/summary_yeast_architecture.txt",row.names = F,quote = F)

df_summary2_lecanorales<- df %>%filter(order=="Lecanorales") %>% group_by(architecture,yeast,presence) %>%   dplyr::summarise(n = n()) %>% dplyr::mutate(freq = 100*n / sum(n)) %>%
  select(-n) %>%filter(architecture!="Unknown") %>% pivot_wider(names_from = presence,values_from=freq )
write.table(df_summary2,"results/tables/summary_yeast_architecture_lecanorales.txt",row.names = F,quote = F)




#make a proportion graph - not polished
df_summary$yeast<-factor(df_summary$yeast, levels = c("cypho", "tremella"), 
                  labels = c("Cyphobasidium", "Tremella"))
df_summary$method<-factor(df_summary$method,levels=c("hit_assembly_abs","hit_reads_abs"),
                          labels=c("In assembly","In reads"))

bars<-ggplot(df_summary,aes(fill=outcome,x=method,y=freq)) +  
  geom_bar(position="fill", stat="identity", colour="black") +
  geom_text(aes(label=scales::percent(df_summary$freq),color=outcome), position=position_fill(vjust=0.5),show.legend = F)+
  scale_y_continuous(name = "", breaks = NULL) +
  scale_fill_manual(values=c("white",col_pres_abs))+
  scale_color_manual(values = c("black","white"))+
  facet_wrap(yeast~.)+
  theme_minimal()+ theme(strip.text.x = element_text(size=12,face = "italic"),axis.title = element_blank())

ggsave("results/figures/metagenome_reanalysis_proportions.png",device="png",height = 10,width=5)


dothist+bars
ggsave("results/figures/lendemer_multi.png",device="png",height = 10,width=15)

#What are non-lecanorales metagenomes with cyphobasidium?
nonlecanorales_cypho_table<-df %>% filter(order!="Lecanorales",yeast=="cypho",presence!="absent") %>% select(c(Lichen.metagenomes,bp,presence,architecture,family, order, class))
write.table(nonlecanorales_cypho_table,"results/tables/nonlecanorales_cypho_table.txt",row.names = F,quote = F, sep="\t")






