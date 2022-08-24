# processing metaxa output
library(tidyverse)
library(stringr)
library(patchwork)
library(pbapply)

matrix_level5<-read.delim("analysis/03_metagenome_reanalysis/metaxa_level_5_combined.txt")
bp<-read.delim("analysis/03_metagenome_reanalysis/bp_report.txt",header=F,col.names=c("metagenome","sequencing_depth"))

#modify the matrix
matrix_level5_long<-matrix_level5 %>% pivot_longer(-Taxa,names_to="sample",values_to="abundance") %>%
  group_by(sample,Taxa) %>% summarize(abundance=sum(abundance)) %>%
  separate(sample,into=c("type","metagenome"),sep="_") %>% mutate(occurrence=ifelse(abundance>0,1,0))


#get most common lineages in assemblies and reads
occ_assembly_level5<-matrix_level5_long %>% filter(type=="assembly")  %>%
  group_by(Taxa) %>% summarize(total_occ=sum(occurrence)) %>% filter(total_occ>0) %>% arrange(desc(total_occ))

occ_reads_level5<-matrix_level5_long %>% filter(type=="reads") %>% mutate(occurrence=ifelse(abundance>0,1,0)) %>%
  group_by(Taxa) %>% summarize(total_occ=sum(occurrence))%>% filter(total_occ>0) %>% arrange(desc(total_occ))

#save the tables
write.table(occ_assembly_level5,"analysis/03_metagenome_reanalysis/occurrence_assembly_metaxa_level5.tsv",sep="\t",quote = F, row.names = F)

write.table(occ_reads_level5,"analysis/03_metagenome_reanalysis/occurrence_reads_metaxa_level5.tsv",sep="\t",quote = F, row.names = F)


#filter
matrix_level5_long_filtered<-matrix_level5_long%>% 
  filter(!grepl('Chloroplast', Taxa)) %>% filter(!grepl('Mitochondria', Taxa)) %>% #removed chloroplasts and mitochondria
  filter(!grepl('Archaea', Taxa)) #removed Archaea sequences which are usually misassigned eukaryotes


#count number of different level5 annotations
lineage_count<-matrix_level5_long_filtered %>% group_by(type,metagenome) %>% filter(abundance>0) %>%
  summarize(n=n()) %>% arrange(desc(n)) %>% pivot_wider(names_from=type,values_from=n) %>%
  left_join(bp)
  
a<-ggplot(lineage_count,aes(x=sequencing_depth,y=assembly)) +
  geom_point()+theme_minimal()+xlab("Sequencing depth")+
  scale_x_continuous(breaks=c(0,10000000000,20000000000,30000000000),labels=c("0","10 Gbp","20 Gbp","30 Gbp"))+
  ylab("Number of lineages \n detected in metagenomic assemblies")+
  geom_smooth(method='gam', formula= y~s(x,bs = "cs"))

r<-ggplot(lineage_count,aes(x=sequencing_depth,y=reads)) +
  geom_point()+theme_minimal()+xlab("Sequencing depth")+
  scale_x_continuous(breaks=c(0,10000000000,20000000000,30000000000),labels=c("0","10 Gbp","20 Gbp","30 Gbp"))+
  ylab("Number of lineages \n detected in metagenomic data")+
  geom_smooth(method='gam', formula= y~s(x,bs = "cs"))



#corresponds with the number of mags
mag_number<-read.delim("analysis/05_MAGs/tables/MAG_counts_summary.tsv")

lineage_count<-lineage_count %>% left_join(mag_number %>%select(metagenome,total))

r_m<-ggplot(lineage_count,aes(x=total,y=reads)) +xlab("Number of MAGs") +
  geom_point()+theme_minimal()+ geom_smooth(method='gam', formula= y~s(x,bs = "cs"))+
  ylab("Number of lineages \n detected in metagenomic data")

a_m<-ggplot(lineage_count,aes(x=total,y=assembly)) +xlab("Number of MAGs") +
  geom_point()+theme_minimal()+ geom_smooth(method='gam', formula= y~s(x,bs = "cs"))+
  ylab("Number of lineages \n detected in metagenomic assemblies")

(a + a_m) / (r + r_m)
ggsave("results/figures/metaxa.png",bg="white",width=10,height=6)

#on avarage, what was the ration of lineages reads/assemblies
lineage_count %>% mutate(read_ass_ratio=reads/assembly) %>% summarize(mean_ratio=mean(read_ass_ratio,na.rm=T))

#how namy lineages in assembly in total
taxa_ass<-matrix_level5_long_filtered %>% select(-metagenome) %>% filter(type=="assembly") %>%  
  filter(abundance>0) %>% select(Taxa) %>% unique()

       

###old

#count total number of rDNA per metagenome
lineage_ab<-matrix_level5_long_filtered %>% group_by(type,metagenome) %>% filter(abundance>0) %>%
  summarize(sum=sum(abundance)) %>% arrange(desc(sum)) %>% pivot_wider(names_from=type,values_from=sum) %>%
  left_join(bp)

a<-ggplot(lineage_ab,aes(x=sequencing_depth,y=assembly)) +
  geom_point()+theme_minimal()

r<-ggplot(lineage_ab,aes(x=sequencing_depth,y=reads)) +
  geom_point()+theme_minimal()
a/r


#remove duplicate: aiming for the number of families
#identify taxonomy that is contained within another entry and remove these rows
#only looked at the assemblies
matrix_level5_long_filtered2<-matrix_level5_long_filtered %>% filter(type=="assembly")

is_contained<-function(x){
  sum(str_detect(as.character(matrix_assembly$Taxa),x))>1
}

matrix_assembly<-matrix_level5_long_filtered2 %>% pivot_wider(names_from=metagenome,values_from=abundance,values_fill=0)

l<-pblapply(as.character(matrix_assembly$Taxa),is_contained)
matrix_assembly$contained<-as.vector(do.call(rbind,l))
matrix_assembly_filtered<-matrix_assembly %>% filter(contained==F) %>% select(-c(occurrence,type,contained))


lineage_number<-matrix_assembly_filtered %>% 
  pivot_longer(-Taxa,values_to="abundance",names_to="metagenome") %>%
  filter(abundance>0) %>% group_by(metagenome) %>% summarize(n=n()) %>%
  arrange(desc(n))





# split into rangs
a <- strsplit(as.character(matrix_level5_long_filtered2$Taxa), ";")
lens <- vapply(a, length, integer(1L)) # or lengths(a) in R 3.2
longdf <- df[rep(seq_along(a), lens),]
longdf$string <- unlist(a)

















is_contained<-function(x){
  sum(str_detect(as.character(matrix_level5_long_filtered$Taxa),x))>1
}

l<-pblapply(as.character(matrix_level5_long_filtered2$Taxa),is_contained)
matrix_level5_long_filtered2$contained<-as.vector(do.call(rbind,l))






#make a heatmap

ggplot(matrix_level5_long%>% filter(type=="reads"), aes(Taxa, metagenome, fill= occurrence)) + 
  geom_tile()




#filter
matrix_level5_2<-matrix_level5%>% distinct() %>% #removed duplicates
  filter(!grepl('Chloroplast', V2)) %>% filter(!grepl('Mitochondria', V2)) %>% #removed chloroplasts and mitochondria
  filter(!grepl('Archaea', V2)) #removed Archaea sequences which are usually misassigned eukaryotes

#identify taxonomy that is contained within another entry and remove these rows
is_contained<-function(x){
sum(str_detect(metaxa$V2,x))>1
}

l<-lapply(metaxa$V2,is_contained)
metaxa$contained<-as.vector(do.call(rbind,l))

metaxa<-metaxa %>% filter(contained==F)

# split into rangs
a <- strsplit(metaxa$V2, ";")
lens <- vapply(a, length, integer(1L)) # or lengths(a) in R 3.2
longdf <- df[rep(seq_along(a), lens),]
longdf$string <- unlist(a)
