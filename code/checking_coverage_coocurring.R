#checking depth of coverage of the two co-occurring algae
library(tidyverse)
library(patchwork)
#setwd("~/Documents/coverage")

mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")

#which metagenomes have both algae?
both_alga<-mags_role %>% select(Genome,metagenome,Lichen.metagenomes) %>% 
  filter(Genome %in% c("public_SRR7232211_concoct_bin.8","private_T1904_concoct_bin.38")) 

both_alga_depth<-mags_role %>% select(Genome,metagenome,depth_cov) %>% 
  filter(Genome %in% c("public_SRR7232211_concoct_bin.8","private_T1904_concoct_bin.38")) %>%
  pivot_wider(names_from=Genome,values_from = depth_cov) %>%
  mutate(cov_ratio=public_SRR7232211_concoct_bin.8/private_T1904_concoct_bin.38) 
  
both_alga_depth$public_SRR7232211_concoct_bin.8[is.na(both_alga_mtg$public_SRR7232211_concoct_bin.8)]<-0
both_alga_depth$private_T1904_concoct_bin.38[is.na(both_alga_mtg$private_T1904_concoct_bin.38)]<-0

write.table(both_alga_depth,"analysis/05_MAGs/tables/two_alga_coverage_depth.tsv",sep="\t",quote = F, row.names = F)

#are they geniunely absent where we think they're absent?
cov_bredth<-read.csv("analysis/05_MAGs/tables/read_mapping/bwa_coverage.csv")
cov_bredth_long<- cov_bredth %>% pivot_longer(-Genome,names_to = "metagenome", values_to="breadth")

both_alga_bredth<-cov_bredth_long %>% filter(Genome %in% colnames(both_alga_depth)) %>%
  pivot_wider(names_from=Genome,values_from = breadth)
##yes they are

#visualize
#breadth of coverage
b<-ggplot(both_alga_bredth,aes(x=public_SRR7232211_concoct_bin.8,y=private_T1904_concoct_bin.38))+
  geom_point() + xlab("public_SRR7232211_concoct_bin.8 \n breadth of coverage")+
  ylab("private_T1904_concoct_bin.38 \n breadth of coverage")+
  theme(plot.margin = unit(c(10,30,10,30), "pt"),
        text = element_text(size=20))

#depth of coverage: only showing the metagenomes that had both

d<-ggplot(both_alga_depth %>% filter(!is.na(cov_ratio)),aes(x=public_SRR7232211_concoct_bin.8,y=private_T1904_concoct_bin.38))+
         geom_point() + geom_abline(intercept = 1)  + xlab("public_SRR7232211_concoct_bin.8 \n depth of coverage")+
  ylab("private_T1904_concoct_bin.38 \n depth of coverage")+
  theme(plot.margin = unit(c(10,30,10,30), "pt"),
        text = element_text(size=20))

b+d
ggsave("results/figures/two_algae_cov.png",bg="white", width=12,height=6,units = "in")



