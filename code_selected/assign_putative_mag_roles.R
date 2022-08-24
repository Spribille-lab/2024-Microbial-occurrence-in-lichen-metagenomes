# Role assignment for MAG/metagenome occurrences
###The idea here that one MAG can play different roles in diff metagenomes, 
###e.g. the same Nostoc can be a phtobiont in one lichen and live in cephalodia in another
###or a lecanoromycete MAG is the mycobiont in one lichen and a low-coverage contaminant in another
###to take into account that, I first make a table here with putative assignments,
###and then manually check this table and cross-reference iit the MAG tree and biology of the source lichen



## 1. Make a table with depths of coverage of each MAG in each metagenome
### this is a ballpark estimation, used the bwa_total table for the number of reads aligned to each MAG and MAG lengths
###only calculated depth for MAGs that are classified as present in a given metagenome (coverage_breadth>50)
library(tidyverse)

###get bwa data
cov_bredth<-read.csv("analysis/05_MAGs/tables/read_mapping/bwa_coverage.csv")
cov_reads<-read.csv("analysis/05_MAGs/tables/read_mapping/bwa_counts-total.csv")

###get mag sizes
mag_size<-read.csv("analysis/05_MAGs/tables/mag_sizes.csv",header=F)
colnames(mag_size)[2]<-"size"
mag_size<-mag_size %>% mutate(Genome = str_replace(V1, ".fa", ""))

###combine all data together
cov_bredth_long<- cov_bredth %>% pivot_longer(-Genome,names_to = "metagenome", values_to="breadth")
cov_reads_long<- cov_reads %>% pivot_longer(-Genome,names_to = "metagenome", values_to="read_count")

data<-cov_bredth_long %>% left_join(cov_reads_long) %>% left_join(mag_size) %>%
  mutate(depth_cov=ifelse(breadth>50,read_count*125/size,0)) %>% select(-V1) %>%
  filter(depth_cov!=0)

###sanity check: are values that are derived from this calculation close to the more robust values given by cmseq?
### We only have cmseq coverage values for MAGs in their "origin" metagenome, but we can use those to check the validity of this approach
cov_check<-read.csv("analysis/05_MAGs/tables/cmseq_coverage.csv") %>% mutate(Genome_tmp = str_replace(bin, ".fa", "")) %>%
  mutate(Genome = str_replace(Genome_tmp, "cov_", "")) %>% select(-Genome_tmp)

data2<-data 
data2$origin_metagenome<-gsub('[a-zA-Z]+_(.*)_.*_.*', '\\1', data2$Genome)
data2<-data2 %>% filter(origin_metagenome==metagenome)

data2<-data2 %>% left_join(cov_check) %>% select(Genome,metagenome,depth_cov, median_cov) %>%
  mutate(diff=abs(median_cov-depth_cov)/depth_cov)
hist(data2$diff)
###most values are within 10% distance from cmseq, the distance is exponentially distributed

## 2. Add taxonomy
mags<-read.delim("analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv") %>%select(-median_cov)
mtg<-read.delim("results/all_metagenome_reanalysis.txt") %>%select(Lichen.metagenomes,Run)

data<-data %>% left_join(mags,by=c("Genome"="mag")) %>% left_join(mtg,by=c("metagenome"="Run"))

data$lineage_broad = "Unknown"
data$lineage_broad[grepl(pattern = "1;131567;2759;33154;4751;", x = data$lineage)] = "Fungi"
data$lineage_broad[grepl(pattern = "1;131567;2759;33090;", x = data$lineage)] = "Chlorophyta"
data$lineage_broad[!is.na(data$bat_bacteria)] = "bacteria_other"
data$lineage_broad[grepl(pattern = "d__Bacteria;p__Cyanobacteria", x = data$bat_bacteria)] = "Cyanobacteria"


## 3. Add putative role assignment
data$putative_role<-"Unknown"
data$putative_role[data$lineage_broad=="Chlorophyta"]<-"photobiont_chloro" #all chlorophyta assumed to be photobionts
data$putative_role[data$lineage_broad=="Cyanobacteria"]<-"photobiont_cyano" #all cyanobacteria assumed to be photobionts, except:

has_alga<-data%>%filter(lineage_broad=="Chlorophyta") %>% select(metagenome) %>% unique() 
data$putative_role[data$lineage_broad=="Cyanobacteria" & data$metagenome %in% has_alga$metagenome]<-"cephalodia_cyano" #except if in a metagenome that has a green alga, in that case they're assumed to be in cephalodia
data$putative_role[data$lineage_broad=="bacteria_other"]<-"bacteria_other"

mtgs_with_many_fungi<-data %>% group_by(metagenome) %>% filter(lineage_broad=="Fungi") %>%
  summarize(n=n(),max_cov = Genome[which.max(depth_cov)]) %>% filter(n>1) #all fungal mags are mycobionts if they are the only fungal mags in a given metagenome
data$putative_role[data$lineage_broad=="Fungi" & !(data$metagenome %in% mtgs_with_many_fungi$metagenome)]<-"mycobiont"
data$putative_role[data$lineage_broad=="Fungi" & data$metagenome %in% mtgs_with_many_fungi$metagenome & data$Genome %in% mtgs_with_many_fungi$max_cov]<-"mycobiont" #otherwise, the highest depth of coverage mag is assigned as mycobiont
data$putative_role[data$lineage_broad=="Fungi" & data$metagenome %in% mtgs_with_many_fungi$metagenome & !(data$Genome %in% mtgs_with_many_fungi$max_cov)]<-"fungi_other" #and the rest are fungi_other


## 4. Save the table
write.table(data,"analysis/05_MAGs/tables/MAG_putative_roles_bwa.tsv",sep="\t",quote = F, row.names = F)


## 5. Analyze the curated table
data_edited<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")

data_edited %>% group_by(metagenome) %>% filter(confirmed_role=="mycobiont_missassigned") %>%
  summarize(n=n(),max_cov = Genome[which.max(depth_cov)]) 
