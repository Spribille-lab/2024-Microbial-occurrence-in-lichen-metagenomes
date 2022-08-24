#building nmds for bacterial occurrence in MAGs

#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(vegan)
library(patchwork)

# 0. Prep
##load metadata on metagenomes
mtg_info<-read.delim("analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt")

##proccess occurrence table
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))
mags_role_filtered<-mags_role %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$confirmed_role=="mycobiont"])


# 1. analyze by MAG, all metagenomes, all bacteria: 
##warning: stress is (nearly) zero: you may have insufficient data
##worked, but all dots except a few outliers are concentrated in one tiny spot
### transform into matrix
bac<-mags_role_filtered%>%filter(lineage_broad %in% c("bacteria_other","Cyanobacteria")) %>%
    select(Genome,bac_genus2,bac_family,bac_order,metagenome)
matrix<-bac %>% group_by(Genome, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=Genome,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)


###plot

data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$Genome <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(Genome,bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

# 2. analyze by genus, all metagenomes, all bacteria: 
##warning: stress is (nearly) zero: you may have insufficient data
##worked, but the plot is weird with 2 metagenomes as complete outliers in both axes
matrix<-bac %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
 # coord_equal() +
  theme_bw()
ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  # coord_equal() +
  theme_bw()

# 3. select only top bacterial lineages (by order) to remove noise, analyze by genus: 
## worked, but the plot is weird with 2 metagenomes as a complete outlier in both axis
bac_top<-bac %>%
  filter(bac_order %in% c("Cyanobacteriales","Rhizobiales","Acetobacterales","Acidobacteriales","Chthoniobacterales","Sphingomonadales"))
  
matrix<-bac_top %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)


###plot

data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  # coord_equal() +
  theme_bw()

# 4. filter metagenomes to remove low-coverage, use all bacteria, analyze by genus:
##failed to converge
depth<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header = F)
colnames(depth)<-c("metagenome","depth")

bac_depth<-bac %>% left_join(depth) %>% filter(depth>2000000000)
# analyze by genus
matrix<-bac_depth %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  # coord_equal() +
  theme_bw()

# 5. filter metagenomes to remove low-coverage, also use only the 5 top bacterial orders, analyze by genus
##failed to converge
bac_top_depth<-bac_depth %>%
  filter(bac_order %in% c("Cyanobacteriales","Rhizobiales","Acetobacterales","Acidobacteriales","Sphingomonadales"))
matrix<-bac_top_depth %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  # coord_equal() +
  theme_bw()

# 6. same as #5 (used top-5 bacterial orders, removed low-coverage), plus remove weird Peltula
##failed to converge

bac_top_depth_remove_pelt<-bac_depth %>%
  filter(bac_order %in% c("Cyanobacteriales","Rhizobiales","Acetobacterales","Acidobacteriales","Sphingomonadales")) %>%
  filter(metagenome!="SRR1531569")

matrix<-bac_top_depth_remove_pelt %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())


ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),alpha=0.5,size=3) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=bac_genus2,colour=bac_family)) +
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

# 7. take top 12 genera and analyze by mag: 
##failed: didn't converge
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())

bac_top_genera<-bac %>%
  filter(bac_genus2 %in% c(annotated_mags$bac_genus2))

matrix<-bac_top_genera %>% group_by(Genome, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=Genome,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=100)

# 8. take top 12 genera and analyze by genus
##failed: didn't converge
matrix<-bac_top_genera %>% group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=3,trymax=200)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=KO),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()

ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),alpha=0.5,size=3) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=bac_genus2))+
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),size=3) + # add the point markers
  #geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2,shape=bac_order,colour=bac_family),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Genome),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()


# 9. Take top families and high-cov metagenomes, analyze by MAG
##failed to converge
matrix<-bac_depth %>% filter(bac_family %in% c("Nostocaceae","Beijerinckiaceae","Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","UBA10450")) %>%
  group_by(Genome, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=Genome,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome

nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=2000)

# 10. Take top families and high-cov metagenomes, analyze by genus
##worked, but Peltula is an outlier very far from other
matrix<-bac_depth %>% filter(bac_family %in% c("Nostocaceae","Beijerinckiaceae","Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","UBA10450")) %>%
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome
nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=1000)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

ggplot() + 
geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),alpha=0.5,size=3) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=bac_genus2))+
   coord_equal() +
  theme_bw()

# 11. Take top families and high-cov metagenomes, analyze by genus,remove Peltula
#finally converged, not great ordinations:
#Stress:     0.1589902 
#Two convergent solutions found after 1435 tries

matrix<-bac_depth %>% filter(metagenome!="SRR1531569",bac_family %in% c("Nostocaceae","Beijerinckiaceae","Acetobacteraceae","Acidobacteriaceae","Sphingomonadaceae","UBA10450")) %>%
  group_by(bac_genus2, metagenome) %>% summarize(n=n()) %>% 
  pivot_wider(names_from=bac_genus2,values_from=n,values_fill=0)
matrix<-matrix %>% data.frame
rownames(matrix)<-matrix$metagenome
nmds=metaMDS(matrix %>% select(-metagenome),k=2,trymax=2000)

###plot
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$metagenome <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores<- data.scores %>% left_join(mtg_info,by=c("metagenome"="Run"))  #  add the grp variable created earlier

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$bac_genus2 <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores %>% left_join(bac %>% select(bac_genus2,bac_family,bac_order) %>% distinct())

sites<-ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=architecture,colour=order),alpha=0.5,size=3) + # add the point markers
 # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=bac_genus2))+
  coord_equal() + xlim(-2,1.3)+ylim(-1,1)+
  theme_bw()
sites_depth<-ggplot() + 
  geom_point(data=data.scores %>% left_join(depth),aes(x=NMDS1,y=NMDS2,shape=architecture,colour=depth),size=3) + # add the point markers
  scale_color_gradient(low="yellow",high="red")+
  coord_equal() + xlim(-2,1.3)+ylim(-1,1)+
  theme_bw()
  
species<-ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=bac_genus2,col=bac_family),position=position_nudge(x = 0.1, y = 0.03))+
  coord_equal() + xlim(-2,1.3)+ylim(-1,1)+
  theme_bw()
sites/species
ggsave("analysis/05_MAGs/exploratory_fig/nmds_sites_species.png",device="png",height = 16,width=14,bg="white")
sites_depth
ggsave("analysis/05_MAGs/exploratory_fig/nmds_sites_depth.png",device="png",height = 9,width=14,bg="white")
