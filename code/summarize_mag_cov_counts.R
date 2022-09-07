# Summarizing MAG coverage and abundance in diff metagenomes


## 1. Get the table (produced by code/assigne_putative_mag_roles.R)
data_edited<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mtg<-read.delim("results/tables/all_metagenome_reanalysis.txt") 


## 2. Analyze total coverage by cetagory
###make a table for each metagenome with combined coverage of diff MAG roles

mag_cov_combined <- data_edited %>% group_by(metagenome,confirmed_role) %>%
  summarise(cov_combined=sum(depth_cov)) %>% 
  pivot_wider(names_from="confirmed_role",values_from = "cov_combined",values_fill = 0) %>%
  mutate(mycobiont_all=mycobiont+mycobiont_missassigned)

###add metadata for the metagenome
mag_cov_combined<-mtg %>% left_join(mag_cov_combined,by=c("Run"="metagenome"))
mag_cov_combined[is.na(mag_cov_combined)]<-0
depth<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F)
colnames(depth)<-c("Run","depth")
mag_cov_combined<-left_join(mag_cov_combined,depth)

### save the table
write.table(mag_cov_combined,"analysis/05_MAGs/tables/MAG_coverage_summary.tsv",sep="\t",quote = F, row.names = F)

#what are the metagenomes that had higher bacterial coverage than mycobiont cov?
high_bacteria<-mag_cov_combined %>% filter(bacteria_other>mycobiont_all & mycobiont_all>0) 
write.table(high_bacteria,"analysis/05_MAGs/tables/MAG_coverage_high_bacteria.tsv",sep="\t",quote = F, row.names = F)

#what are the metagenomes that had higher algal coverage than mycobiont cov?
high_alga<-mag_cov_combined %>% filter(photobiont_chloro>mycobiont_all & mycobiont_all>0) 
write.table(high_alga,"analysis/05_MAGs/tables/MAG_coverage_high_alga.tsv",sep="\t",quote = F, row.names = F)


## 3. Count MAGs in metagenomes
mag_number<- data_edited %>% group_by(metagenome,confirmed_role) %>%
  summarise(mags_in_a_group=n()) %>% 
  pivot_wider(names_from="confirmed_role",values_from = "mags_in_a_group",values_fill = 0) %>%
  mutate(total=bacteria_other+ cypho + fungi_other + mycobiont+mycobiont_missassigned + photobiont_chloro + photobiont_cyano + tremella+cephalodia_cyano+algae_other) 
###add metagenomes that didn't produce any MAG
no_mag_mtg<-data.frame(mtg$Run[!(mtg$Run %in% mag_number$metagenome)])
colnames(no_mag_mtg)<-"metagenome"
mag_number<-mag_number %>% plyr:::rbind.fill(no_mag_mtg)
mag_number[is.na(mag_number)] <- 0

##add extra info
mag_number <- mag_number %>% left_join(mtg,by=c("metagenome"="Run")) %>% left_join(depth,by=c("metagenome"="Run")) 

write.table(mag_number,"results/tables/MAG_counts_summary.tsv",sep="\t",quote = F, row.names = F)

### Visualize MAG numbers

ggplot(mag_number,aes(x=depth,y=total))+geom_point(aes(color=architecture),alpha=0.5) + 
  theme_minimal()+ theme(aspect.ratio=1)+ geom_smooth(method='gam', formula= y~s(x,bs = "cs"))+
  ylab("Number of MAGs")+xlab("Sequencing depth")+
  scale_x_continuous(breaks=c(0,10000000000,20000000000,30000000000),labels=c("0","10 Gbp","20 Gbp","30 Gbp"))+
ggsave("results/figures/mags_vs_depth.png",bg="white")

ggplot(mag_number)+geom_point(aes(x=depth,y=total,color=architecture),alpha=0.5) + 
  theme_minimal()+ theme(aspect.ratio=1)+scale_x_log10()
ggsave("analysis/05_MAGs/exploratory_fig/mags_vs_depth_log.png",bg="white")


### fitting a curve
lin<-lm(total~depth,mag_number)
log<-lm(total~log(depth),mag_number)
