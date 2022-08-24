#prepare fastas for dNdS analysis
#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(Biostrings)
library(stringr)
library(ips)

# 1. Select KOs to check
calvin_kos<-c("K00855", "K01601","K01602", "K00927", "K00150","K00134","K01623","K01624",
    "K03841","K02446","K11532","K00615","K01623","K01624","K11532","K00615","K01807","K01808")

photosystem_kos<-c("K08928","K08929")

bchlorophyll_kos<-c("K04035","K19073","K04037","K04038","K04039","K11336","K11334",
    "K11335","K11333","K11337","K04040","K10960")

ppp_kos<-c("K13937","K00036","K19243","K01057","K07404","K00033", "K01783", "K01807",'K01808',
  "K00615", "K00616", "K01810","K06859","K13810",'K15916')
tonb<-c("K03832","K03561","K03559")

# 2. Add info about mags
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())

# 3. Combine with KO table
kegg_combined<-read.delim("analysis/07_annotate_MAGs/summarized_outputs/kegg_combined.txt",col.names = c("locus","KO","Genome"))

kegg_combined<-kegg_combined %>% left_join(annotated_mags) %>% 
  mutate(Pathway=ifelse(KO %in% calvin_kos,"Calvin",ifelse(KO %in% photosystem_kos,"Photosystem", ifelse(KO %in% bchlorophyll_kos, "bchlorophyll",ifelse(KO %in% ppp_kos, "PPP",ifelse(KO %in% tonb, "tonb",NA))))))

selected_loci<-kegg_combined %>% filter(!is.na(Pathway),bac_family %in% c("Acetobacteraceae","Beijerinckiaceae"))

# 4. Read fastas
read_fasta<-function(mag_list,suffix){
  
  read_fasta_by_mag<-function(mag,suffix){
    filename<-paste0("analysis/07_annotate_MAGs/",mag,"/",mag,suffix)
    if (suffix==".faa"){
      fasta<-readAAStringSet(filename)
    }else{
      fasta<-readDNAStringSet(filename)
    }
    return(fasta)}
  
  l<-lapply(mag_list,read_fasta_by_mag,suffix=suffix)
  fasta_all<-do.call(c,l)
  return(fasta_all)   
}

all_faa<-read_fasta(unique(selected_loci$Genome),".faa")
names(all_faa)<-names(all_faa)%>% word(1)

all_ffn<-read_fasta(unique(selected_loci$Genome),".ffn")
names(all_ffn)<-names(all_ffn)%>% word(1)


# 5. Initial dnds
dnds_by_ko_taxa<-function(ko,taxa,taxa_level){
  taxa_level<-enquo(taxa_level)
  loci<-selected_loci %>% filter(KO==ko,!!taxa_level==taxa)
  ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
  ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
  ffn_aligned<-ffn_aligned %>% del.colgapsonly(threshold = 0.00000001)
  dnds_matrix<-dnds(ffn_aligned) 
  df<-summary(na.omit(dnds_matrix)) %>% as.vector() %>% t() %>% data.frame()
  colnames(df)<-c("Min","1st.Qu","Median","Mean", "3rd.Qu","Max")
  df2<-data.frame("KO"=ko,"taxa"=taxa)
  df3<-cbind(df2,df)
  return(df3)
}


ko_list<-tonb

loci<-selected_loci %>% filter(KO %in% ko_list,bac_genus2=="Acetobacteraceae gen. sp.")
loci2<-loci %>% group_by(KO) %>% summarize(n=n()) %>% filter(n>1)
ko_list2<-loci2$KO

l<-lapply(ko_list2,dnds_by_ko_taxa,"Acetobacteraceae gen. sp.",bac_genus2)
do.call(rbind,l)


# 6. Follow up on KOs with dNdS > 1. 
#6a
ko<-"K04040"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="CAHJXG01")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #identified gaps
ffn_trimmed<-ffn_aligned[,55:929]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

#6b
ko<-"K04040"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Acetobacteraceae gen. sp.")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #identified gaps
ffn_trimmed<-ffn_aligned[,224:979]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed,details = T) 
summary(dnds_matrix)
# dNdS < 1 for one pair that weren't NaN -> probably exonerated, will look more into

#6c
ko<-"K11333"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="CAHJXG01")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #identified gaps
ffn_trimmed<-ffn_aligned[,100:1017]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

#6d
ko<-"K04038"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Lichenihabitans")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #identified gaps
ffn_trimmed<-ffn_aligned[,46:1308]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

#6e
ko<-"K04038"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="RH-AL1")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #identified gaps
ffn_trimmed<-ffn_aligned[,602:1324]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

#6f
ko<-"K00615"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="RH-AL1")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated

loci1<-loci[c(1,3,4),]
ffn<-all_ffn[names(all_ffn) %in% loci1$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #no gaps
dnds_matrix<-dnds(ffn_aligned) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

loci2<-loci[c(2,5),]
ffn<-all_ffn[names(all_ffn) %in% loci2$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,440:700]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

#6g
ko<-"K00615"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Lichenihabitans")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated

loci1<-loci[c(1,3:6,8),]
ffn<-all_ffn[names(all_ffn) %in% loci1$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,550:999]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

loci2<-loci[c(2,7),]
ffn<-all_ffn[names(all_ffn) %in% loci2$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,199:651]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated


#6h
ko<-"K00615"

loci<-selected_loci %>% filter(KO==ko,bac_genus2=="CAHJXG01")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated

loci1<-loci[c(1:2,4:7,9:12),]
ffn<-all_ffn[names(all_ffn) %in% loci1$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,1451:1750]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed,details=T) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

loci2<-loci[c(3,8),]
ffn<-all_ffn[names(all_ffn) %in% loci2$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,52:804]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

# 6i
ko<-"K03559"
loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Acetobacteraceae gen. sp.")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated
ffn_trimmed<-ffn_aligned[,1685:1793]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
dnds_matrix<-dnds(ffn_trimmed,details=T) 

# 6j
ko<-"K03561"
loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Acetobacteraceae gen. sp.")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated

loci1<-loci[c(1,4,5,7),]
ffn<-all_ffn[names(all_ffn) %in% loci1$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,31:660]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated

# 6k 
ko<-"K03832"
loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Acetobacteraceae gen. sp.")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated
ffn_trimmed<-ffn_aligned[,1100:1150]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS  = 0.96!!!

# 6l 
ko<-"K03832"
loci<-selected_loci %>% filter(KO==ko,bac_genus2=="CAHJXG01")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated

# 6m 
ko<-"K03832"
loci<-selected_loci %>% filter(KO==ko,bac_genus2=="Lichenihabitans")
ffn<-all_ffn[names(all_ffn) %in% loci$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #look like 2 orthologs need to be separated
ffn_trimmed<-ffn_aligned[,448:550]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)

loci1<-loci[c(1,3),]
ffn<-all_ffn[names(all_ffn) %in% loci1$locus ]
ffn_aligned<-mafft(as.DNAbin(ffn),maxiterate=1000,method="genafpair")
checkAlignment(ffn_aligned) #gaps
ffn_trimmed<-ffn_aligned[,490:590]
checkAlignment(ffn_trimmed)
dnds_matrix<-dnds(ffn_trimmed) 
summary(dnds_matrix)
# dNdS <<1 -> exonerated


