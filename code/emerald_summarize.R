##summarizing the emeraldBGC results
library(dplyr)
library(tidyr)
source("code/utils.R")


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
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")


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



## save human-readable table
long_table<-emerald_ann %>% select(sample, bac_genus2,bac_family,bac_order,cluster_id, mibig_class, jaccard,nearest_mibig)
write.table( long_table, "results/tables/bgc_all_good_hits.txt", sep='\t',quote = F, row.names = F, col.names = T)
