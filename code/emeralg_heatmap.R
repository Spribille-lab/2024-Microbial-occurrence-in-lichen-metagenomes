library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
source("code/utils.R")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(RColorBrewer)
options(repr.plot.width=20, repr.plot.height=40)
theme_set(theme_minimal(base_size = 23))


#defining functions

files_from_folder <- function(dirpath){
  filenames <- list.files(path = dirpath)
  listlength <- length(filenames)
  seq_list <- vector("list", listlength)
  for (i in 1:listlength){
    sample <- filenames[i]
    seq_path <- paste(dirpath, sample, "/", sample, ".gbk.emerald/", sample, ".gbk.emerald.full.gff",sep = "")
    if (file.exists(seq_path)) {
      if (file.size(seq_path) > 16) {
        seq_input <- read.gff(seq_path)
      } else {
        seq_input <- data.frame(seqid = "N/A", source = "N/A", type = "N/A", start = "N/A", end = "N/A", score = "N/A", strand = "N/A", phase = "N/A", attributes = "N/A")
      }
    }
    seq_list[[i]] <- seq_input
  }
  
  names(seq_list) <- filenames
  for (i in seq_along(seq_list)) seq_list[[i]]$sample <- filenames[i]
  bgc <- Reduce(rbind, seq_list)
  
  bgc <- bgc %>% separate(attributes, c("cluster_id", "nearest_mibig", "mibig_class", "jaccard", "partial_cluster"), sep = ";")
  
  for (i in 1:nrow(bgc)){
    for (j in 9:13){
      bgc[i,j] <- gsub(".*=", "", bgc[i,j])
    }
  }
  bgc
}

emeraldbgc <- files_from_folder("analysis/07_annotate_MAGs/emeraldbgc/")


###get info on mags
## 2. read data
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
annotated_mags<-annotated_mags %>% filter(Genome!="public_SRR14722130_metawrap_bin.2" & Genome!="private_T1889_metawrap_bin.7")


## 3. get taxonomy
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")
mags_role$bac_class <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="c")

mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

emeraldbgc <-emeraldbgc %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order,bac_class,bac_phylum) %>% distinct(), by=c("sample"="Genome"))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())
annotated_mags_arranged<-annotated_mags %>% 
  group_by(bac_family) %>%
  arrange(bac_genus2, .by_group=T)

##subset of only mags of interest. remove hits with jaccard values >0.7
emerald_ann<-emeraldbgc %>% filter(sample %in% annotated_mags$Genome, jaccard<0.7)

## look at individual classes
terp<-emerald_ann %>% filter(mibig_class=="Terpene") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

polyket<-emerald_ann %>% filter(mibig_class=="Polyketide") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

nrp<-emerald_ann %>% filter(mibig_class=="NRP") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

ripp<-emerald_ann %>% filter(mibig_class=="RiPP") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

sac<-emerald_ann %>% filter(mibig_class=="Saccharide") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

other<-emerald_ann %>% filter(mibig_class=="Other") %>% group_by(bac_genus2,nearest_mibig) %>%
  summarize(n=n()) %>% pivot_wider(values_from = n,values_fill=0,names_from = bac_genus2)

### visualize BGCs of interest
bgc_of_interest<-read.delim('analysis/07_annotate_MAGs/bgc_of_interest.txt')
bgc_of_interest$new_name<-paste0(bgc_of_interest$mibig,": ",bgc_of_interest$compaund," (",bgc_of_interest$function.,")")
bgc_of_interest <- bgc_of_interest %>% filter(mibig %in%emerald_ann$nearest_mibig )

df<- emerald_ann  %>%
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill=0) %>%
  select(sample,one_of(bgc_of_interest$mibig))


##prepare matrix
M <- t(as.matrix(df[,2:ncol(df)]))
colnames(M) = df$sample
new_rownames<-data.frame(rownames(M)) %>% left_join(bgc_of_interest,by=c("rownames.M."="mibig"))
rownames(M)<-new_rownames$new_name

###set colors
genus_node_color<-c('Acetobacteraceae gen. sp.' = brewer.pal(12,"Set3")[2],
                    "CAHJXG01"  = brewer.pal(12,"Set3")[3],
                    "CAIMSN01" = brewer.pal(12,"Set3")[4],
                    "LMUY01"  = brewer.pal(12,"Set3")[5],
                    "VCDI01" = brewer.pal(12,"Set3")[6],
                    "CAHJWL01"  = brewer.pal(12,"Set3")[7],
                    "EB88"  = brewer.pal(12,"Set3")[8],
                    "Terriglobus" = brewer.pal(12,"Set3")[9],
                    "Lichenihabitans" = brewer.pal(12,"Set3")[10],
                    "RH-AL1" = brewer.pal(12,"Set3")[11],
                    "Nostoc" = brewer.pal(12,"Set3")[1],
                    "Sphingomonas" = brewer.pal(12,"Set3")[12],
                    "CAHJWO01" = brewer.pal(8,"Dark2")[1])
family_node_color<-c("Acetobacteraceae" = brewer.pal(8,"Dark2")[2],
                     "Acidobacteriaceae" = brewer.pal(8,"Dark2")[3],
                     "Beijerinckiaceae"  = brewer.pal(8,"Dark2")[4],
                     "Nostocaceae" = brewer.pal(8,"Dark2")[5],
                     "Sphingomonadaceae"  = brewer.pal(8,"Dark2")[6],
                     "UBA10450" = brewer.pal(8,"Dark2")[7])

#make top annotation for MAGs
genus_df<-data.frame(colnames(M)) %>% left_join(annotated_mags %>%select(Genome, bac_genus2) %>%unique(),by=c("colnames.M."="Genome"))
genus<-genus_df$bac_genus2

family_df<-data.frame(colnames(M)) %>% left_join(annotated_mags %>%select(Genome, bac_family) %>%unique(),by=c("colnames.M."="Genome"))
family<-family_df$bac_family

genus_ordered<-annotated_mags_arranged$bac_genus2 %>% unique()


ta = HeatmapAnnotation(df = data.frame(family=family,genus = genus), show_annotation_name=T,
                       annotation_legend_param = list(genus=list(labels=genus_ordered,at=genus_ordered),direction = "vertical"),
                       col = list(genus = genus_node_color,family=family_node_color)
)


###row annotation
class <- new_rownames %>% left_join(bgc_of_interest %>% ungroup(),by=c("rownames.M."="mibig"))
class <-class$function..y
class_ordered<- bgc_of_interest$function..y %>% unique()

class_color<-c( "Carotenoid"  = brewer.pal(8,"Pastel2")[1],
               "Antibiotic"  = brewer.pal(8,"Pastel2")[2],
               "Cytotoxic"  = brewer.pal(8,"Pastel2")[3],
               "Hormone"  = brewer.pal(8,"Pastel2")[4],
               "Siderophore"  = brewer.pal(8,"Pastel2")[5],
               "Exopolysaccharide"  = brewer.pal(8,"Pastel2")[6],
               "Biosurfactant"  = brewer.pal(8,"Pastel2")[7],
               "Phytotoxic"  = brewer.pal(8,"Pastel2")[8],
               "Sunscreen"  = brewer.pal(8,"Dark2")[1],
               "Antibiotic Cytotoxic"  = brewer.pal(8,"Dark2")[2],
               "Antimicrobial"  = brewer.pal(8,"Dark2")[3],
               "Antifungal"  = brewer.pal(8,"Dark2")[4],
               "Cell Wall Polysacchharide"  = brewer.pal(8,"Dark2")[5]
               )

ra = rowAnnotation(df=data.frame(class  = class ),
                   col = list(class  = class_color),
                   annotation_legend_param = list(class =list(labels=class_ordered,at=class_ordered,title = "function")))


### set color scheme
tmp<-df %>% pivot_longer(-sample,names_to="family",values_to="count")
col_fun = colorRamp2(c(0,1, max(tmp$count)), c("white","orange", "red"))


hm<-Heatmap(M, show_column_names = F, name = " ",
            row_order = bgc_of_interest$new_name[bgc_of_interest$new_name %in% rownames(M)],
            column_order=annotated_mags_arranged$Genome[annotated_mags_arranged$Genome %in% colnames(M)],
            top_annotation = ta,
            right_annotation = ra,
            col = col_fun,
            row_names_side = "left",
            heatmap_legend_param = list(
              title = "# of Genes", at = c(0, 1, 15, 30, 45), 
              labels = c(0, 1, 15, 30, 45))
)
hm

pdf(file="analysis/07_annotate_MAGs/summarized_outputs/bgc_heatmap.pdf",width=17,height=30)
draw(hm,annotation_legend_side = "top")
dev.off()

## save human-readable table
hr_table<-emerald_ann %>% group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0)  %>%
   left_join(annotated_mags,by=c("sample"="Genome")) 

hr_table<-hr_table %>%  select(sort(current_vars())) %>% 
  relocate(sample,bac_genus2,bac_family,bac_order) %>% group_by(bac_family) %>%
  arrange(bac_genus2, .by_group=T)
write.table( hr_table, "analysis/07_annotate_MAGs/summarized_outputs/bgc_summarized.txt", sep='\t',quote = F, row.names = F, col.names = T)


##how many mags have exopol bgs?
eps_bgc_list<-bgc_of_interest %>% filter(function.=="Exopolysaccharide")

#how many eps clusters per mag
eps_by_mag<-emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,all_of(eps_bgc_list$mibig)) %>%
  pivot_longer(-sample,names_to="nearest_mibig",values_to="n") %>%
  group_by(sample) %>% summarise(n_sacc=sum(n)) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome"))

#average by family
eps_by_fam<-emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,all_of(eps_bgc_list$mibig)) %>%
  pivot_longer(-sample,names_to="nearest_mibig",values_to="n") %>%
  group_by(sample) %>% summarise(n_sacc=sum(n)) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome")) %>%
  group_by(bac_family) %>%  summarize(avg=mean(n_sacc))

# % of mags with eps clusters (out of the 63 annotated)
emerald_ann %>% filter(nearest_mibig %in% eps_bgc_list$mibig) %>%
  select(sample) %>% distinct() %>% nrow()
#  % of mags with eps clusters (out of the 63 annotated), grouped by genus
emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,all_of(eps_bgc_list$mibig)) %>%
  pivot_longer(-sample,names_to="nearest_mibig",values_to="n") %>%
  group_by(sample) %>% summarise(n_sacc=sum(n)) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome")) %>%
  mutate(presence=ifelse(n_sacc>0,1,0)) %>%
  group_by(bac_genus2,presence) %>% summarise(n=n()) %>% 
  pivot_wider(values_from = n,names_from = presence,values_fill = 0)

#  % of mags with eps clusters (out of the 63 annotated), grouped by family
emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,all_of(eps_bgc_list$mibig)) %>%
  pivot_longer(-sample,names_to="nearest_mibig",values_to="n") %>%
  group_by(sample) %>% summarise(n_sacc=sum(n)) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome")) %>%
  mutate(presence=ifelse(n_sacc>0,1,0)) %>%
  group_by(bac_family,presence) %>% summarise(n=n()) %>% 
  pivot_wider(values_from = n,names_from = presence,values_fill = 0)







##count saccharide bgcs pre genome
sacch_mibig<-emerald_ann %>% filter(mibig_class=="Saccharide") %>% select(nearest_mibig) %>% distinct()
t<-emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,all_of(sacch_mibig$nearest_mibig)) %>%
  pivot_longer(-sample,names_to="nearest_mibig",values_to="n") %>%
  group_by(sample) %>% summarise(n_sacc=sum(n)) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome"))
  


emerald_ann %>% filter(mibig_class=="Saccharide") %>% group_by(sample) %>%
  summarize(n=n()) %>% 
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome")) %>%
  ###average by family
  group_by(bac_family) %>%  summarize(avg=mean(n))

##which mags have BGC0000731
emerald_ann %>% 
  group_by(nearest_mibig,sample) %>% summarize(n=n()) %>%
  pivot_wider(names_from=nearest_mibig,values_from=n,values_fill = 0) %>% 
  select(sample,BGC0000731) %>%
  left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct(),by=c("sample"="Genome")) %>%
  group_by(bac_family) %>%  summarize(avg=mean(BGC0000731))
  