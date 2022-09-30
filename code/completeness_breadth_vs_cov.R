# mag completenes as a function of coverage
library(tidyverse)

#load QC results

df<-read.delim("analysis/05_MAGs/tables/euks_cov_QC.txt")
euk_compl<-ggplot(df,aes(x=median_cov,y=completeness))+geom_point()+
  scale_x_log10()+theme_minimal()+xlab("Median depth of coverage")+
  ylab("Completeness according to EukCC, %")


#what was the lowest coverage as estimated by bwa that resuted in a highly complete MAG?
df_filtered<-df %>% filter(completeness>95) #only include MAGs with completenes high according to EukCC
df_filtered$Genome<-str_replace(df_filtered$genome,".fa","")
  
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mags_role_filtered<-mags_role %>% filter(breadth>95)

mags_role_filtered %>% inner_join(df_filtered) %>% 
  select(Genome, metagenome, completeness,depth_cov,breadth) %>% 
  arrange(depth_cov) #the lowest depth_cov is 4

#is it normal for a 4X genome to have a good breath of coverage?
low_cov<-mags_role %>% inner_join(df_filtered) %>% 
  select(Genome, metagenome, completeness,depth_cov,breadth) %>%
  filter(depth_cov<10 )

breadth<-ggplot(low_cov,aes(x=depth_cov,y=breadth))+geom_point()+
  theme_minimal()+xlab("Median depth of coverage")+
  ylab("Breadth of coverage, %")

euk_compl + breadth
ggsave("results/figures/completenes_and_breadth_vs_cov.png",bg="white",width=10,height=5)


