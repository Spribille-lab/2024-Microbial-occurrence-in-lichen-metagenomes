source("code/utils.R")
library(tidyverse)
library(patchwork)

##load the mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

##add labels
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$label<-mags_role$bac_family
mags_role$label[mags_role$confirmed_role=="mycobiont"]<-"Mycobiont"
mags_role$label[mags_role$confirmed_role=="photobiont_chloro"]<-"Trebouxiophyceae"
mags_role$label[mags_role$confirmed_role=="cypho"]<-"Cystobasidiomycetes"
mags_role$label[mags_role$confirmed_role=="tremella"]<-"Tremellomycetes"

##define the list of metagenomes that yeilded at least one mag
metagenomes_with_mags<-mags_role$metagenome %>% unique


##make a table with mag occurrences of the groups of interest
groups_to_include<-c("Trebouxiophyceae","Acidobacteriaceae",
                     "Beijerinckiaceae","Acetobacteraceae",
                     "Sphingomonadaceae","Cystobasidiomycetes","Tremellomycetes")

mag_presence <- mags_role %>% dplyr::select(metagenome,breadth,label) %>%
  dplyr::filter(label %in% groups_to_include) %>% 
  dplyr::filter(breadth>50) %>% dplyr::select(-breadth) %>% distinct() %>%
  group_by(metagenome,label) %>% summarize(n=n()) %>%
  pivot_wider(values_from = n, values_fill = 0, names_from = label) 

##add empty lines for the metagenomes that didn't have any mags of interest (but still had at least one mag)
mag_presence <- data_frame("metagenome"=metagenomes_with_mags) %>% left_join(mag_presence) %>%
  mutate_if(is.numeric, ~replace_na(.,0)) %>%
  pivot_longer(-metagenome,names_to="group",values_to="presence_mag")

##add depth info
depth<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F)
colnames(depth)<-c("metagenome","depth")
mag_presence <- mag_presence %>% left_join(depth)

##add results of rRNA screeining
rdna<-read.delim2("analysis/03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv")
rdna <- rdna %>% filter(type=="reads") %>% select(-type) %>% pivot_longer(-metagenome,names_to="group",values_to="presence_reads")
  
  

#plot
df<-mag_presence %>% left_join(rdna)%>%
  mutate(presence=ifelse(presence_mag>0,"MAG detected",
                         ifelse(presence_reads>0, "rRNA detected","not detected")),
         presence_dummy=ifelse(presence_mag>0,0,
                         ifelse(presence_reads>0, 1,2))) %>%
  arrange(across(.cols=c("presence_dummy","depth"))) %>% 
  rowid_to_column()
df$presence<-factor(df$presence,levels = c("MAG detected","rRNA detected","not detected"))

p1<-ggplot(df %>% filter(group=="Beijerinckiaceae"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+ guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p2<-ggplot(df %>% filter(group=="Acetobacteraceae"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+ guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p3<-ggplot(df %>% filter(group=="Acidobacteriaceae"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+ guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p4<-ggplot(df %>% filter(group=="Sphingomonadaceae"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("sequencing depth, bp")+ guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p5<-ggplot(df %>% filter(group=="Trebouxiophyceae"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+ guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p6<-ggplot(df %>% filter(group=="Tremellomycetes"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+guides(fill="none")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())

p7<-ggplot(df %>% filter(group=="Cystobasidiomycetes"))  +
  geom_bar(aes(x =  reorder(metagenome, rowid),  y = depth, fill = presence), stat = "identity",  width = 0.8)+
  ylab("")+
  scale_fill_manual(values = c("#28620f","#5dbd33","#c9edba"))+
  theme(legend.position="bottom",
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())

p5 / p2 / p1/ p3 / p4 / p6 / p7
ggsave("results/taxa_presence_cov.svg",width=120,height=115,device="svg",bg="white",units="mm")

#save as a table
df2 <- df %>% select(metagenome,group,depth,presence)
write.table(df2,"results/tables/group_presence_depth.tsv",sep="\t",quote = F, row.names = F,col.names = T)


