#setwd("~/Documents/coverage")
library(tidyverse)
library(ape)
library(phytools)
library(stringr)
library(RColorBrewer)
source("code/utils.R")

## 1. Combine all info about Rhizobiales MAGs in one place

### GenBank data
genbank_table<-read.csv("analysis/09_rhizobiales_phylogeny/Table_1.csv",header=T,skip=1,
                        na.strings="")

genbank_table<-genbank_table %>% filter(!is.na(assembly_accession)) %>%
  dplyr::select(assembly_accession,organism_name)

### GTDBTk classify results
gtdbtk<-read.delim("analysis/09_rhizobiales_phylogeny/gtdbtk_classify/gtdbtk.bac120.summary.tsv",header=T)

gtdbtk<-gtdbtk %>% mutate(source= ifelse(grepl("GCF",user_genome),"genbank","own")) %>%
  mutate(Genome = ifelse(grepl("GCF",user_genome),gsub("([^_]*_[^_]*)_.*$", "\\1", user_genome),str_replace(user_genome,".fa",""))) %>%
  mutate(bac_order=sapply(pplacer_taxonomy, gtdb_get_clade, clade="o")) %>%
  mutate(bac_family=sapply(pplacer_taxonomy, gtdb_get_clade, clade="f")) %>%
  mutate(bac_genus=sapply(pplacer_taxonomy, gtdb_get_clade, clade="g")) %>%
  mutate(old_name=user_genome)%>%
  dplyr::select(Genome,bac_order,bac_family,bac_genus,source,old_name)

#### remove _A from order/family names
gtdbtk$bac_order2<-gsub("([^_]*)_.*$", "\\1",gtdbtk$bac_order)
gtdbtk$bac_family2<-gsub("([^_]*)_.*$", "\\1",gtdbtk$bac_family)

#### add info from the table from paper (with names of these bacteria)
df<-gtdbtk %>% left_join(genbank_table,by=c("Genome"="assembly_accession")) %>%
  mutate(new_label = ifelse(source=="genbank",paste(Genome,organism_name,sep=" "),paste(Genome,bac_genus,sep=" ")))




## 2. Searching for NifH genes

### results of grep serch of "nitrogenase iron protein" in the ncbi annotation of genomes
grep<-read.delim("analysis/09_rhizobiales_phylogeny/grep_nitrogenase_iron_protein.txt",header=F)
grep$assembly_accession<-gsub("^.*/([^_]*_[^_]*)_.*$", "\\1", grep$V1)

list_grep_nifh<-unique(grep$assembly_accession)

### results of blast search
blast<-read.delim("analysis/09_rhizobiales_phylogeny/blast_nifh/blast_nifh.txt",header=F)
colnames(blast)<-c("file","method","results")
blast$assembly_accession<-gsub("([^_]*_[^_]*)_.*$", "\\1", blast$file)


### are tblastn and blastp results consistent?
### only have blastp search for NCBI genomes

blast_ncbi<-blast %>% filter(grepl("GCF",assembly_accession))
donot_have_both<-blast_ncbi %>% group_by(assembly_accession) %>% summarize(n=n()) %>% filter(n!=2)

blast_summary<-blast_ncbi %>% select(-file) %>% 
  filter(!(assembly_accession %in% donot_have_both$assembly_accession)) %>%
  pivot_wider(values_from = results, names_from=method) %>%
  mutate(is_consistent = ifelse(tblastn==blastp,T,F))

blast_summary %>% filter(is_consistent == F)
 ####result: is consistent! no difference between tblastn and blastp results

### are blast results and grep results based on the ncbi annotations consistent?
blast_summary<-blast_summary %>% mutate(grep_results = ifelse(assembly_accession %in% list_grep_nifh,1,0)) %>%
  mutate(grep_consistent = ifelse(grep_results==tblastn,T,F)) 

blast_summary %>% filter(grep_consistent == F)

 ####one inconsistency: GCF_002879535.1 has it according to blast, but not according to grep. in annotations this protein si nitrogenase reductase, partial 

### prepared table for presence/absence of NifH based on tblastn search

nifh_presence<-blast %>% filter(method=="tblastn") %>%
  mutate(Genome = ifelse(grepl("GCF",assembly_accession),assembly_accession,str_replace(file,".fa",""))) %>%
  select(Genome,results)

df<-df %>% left_join(nifh_presence)


## 3. search for methane metabolism genes
####read blast results
blast_pmoc<-read.delim("analysis/09_rhizobiales_phylogeny/blast_methane/blast_methane.txt",header=F)
colnames(blast_pmoc)<-c("file","method","presence")
blast_pmoc$gene<-"pmoc"

blast_mmox<-read.delim("analysis/09_rhizobiales_phylogeny/blast_methane/blast_mmox.txt",header=F)
colnames(blast_mmox)<-c("file","method","presence")
blast_mmox$gene<-"mmox"

blast_gmas<-read.delim("analysis/09_rhizobiales_phylogeny/blast_methane/blast_gmas.txt",header=F)
colnames(blast_gmas)<-c("file","method","presence")
blast_gmas$gene<-"gmas"

blast_xoxf<-read.delim("analysis/09_rhizobiales_phylogeny/blast_methane/blast_xoxf.txt",header=F)
colnames(blast_xoxf)<-c("file","method","presence")
blast_xoxf$gene<-"xoxf"

blast_mxaf<-read.delim("analysis/09_rhizobiales_phylogeny/blast_methane/blast_mxaf.txt",header=F)
colnames(blast_mxaf)<-c("file","method","presence")
blast_mxaf$gene<-"mxaf"

blast_methane<-rbind(blast_mxaf,blast_xoxf,blast_gmas,blast_mmox,blast_pmoc)

##### process the file
blast_methane$assembly_accession<-gsub("([^_]*_[^_]*)_.*$", "\\1", blast_methane$file)

methane_presence<-blast_methane %>% 
  mutate(Genome = ifelse(grepl("GCF",assembly_accession),assembly_accession,str_replace(file,".fa",""))) %>%
  dplyr::select(Genome,presence,gene) %>%
  pivot_wider(values_from=presence,names_from=gene)
  
df<-df %>% left_join(methane_presence)


## 4. read and rename tree
# read bacterial tree
tree<-read.newick("analysis/09_rhizobiales_phylogeny/iqtree/rhizobiales.contree")
old_labels<-data.frame(tree$tip.label)


#make new labels
lables_both<-left_join(old_labels,df,by=c("tree.tip.label"="old_name"))
tree$tip.label<-lables_both$new_label
write.tree(tree, file='analysis/09_rhizobiales_phylogeny/iqtree/rhizobiales.contree.renamed')



## 5. make annotation files

### family
family_df<-lables_both %>% select(new_label,bac_family2)

family_summary <- family_df %>% group_by(bac_family2) %>% summarize(n=n())

family_df<-family_df %>% mutate(col=ifelse(bac_family2=="Rhizobiaceae",brewer.pal(8,"Pastel2")[1],
        ifelse(bac_family2=="Beijerinckiaceae",brewer.pal(8,"Pastel2")[2],
               ifelse(bac_family2=="Xanthobacteraceae",brewer.pal(8,"Pastel2")[3],
                      ifelse(bac_family2=="Devosiaceae",brewer.pal(8,"Pastel2")[4],
                             ifelse(bac_family2=="Hyphomicrobiaceae",brewer.pal(8,"Pastel2")[5],
                                    ifelse(bac_family2=="Aestuariivirgaceae",brewer.pal(8,"Pastel2")[6],
                                           ifelse(bac_family2=="Pleomorphomonadaceae",brewer.pal(8,"Pastel2")[7],
                                                  ifelse(bac_family2=="Methyloligellaceae",brewer.pal(8,"Pastel2")[8],"") ))))))))
family_df$new_label<-str_replace_all(family_df$new_label," ","_")

itol_fam<-family_df[,c(1,3,2)]


cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,Family\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_fam.txt")
write.table(itol_fam,"analysis/09_rhizobiales_phylogeny/iqtree/itol_fam.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)

itol_fam2<-itol_fam %>% mutate(type="range")
itol_fam2<-itol_fam2[,c(1,4,2,3)]
cat("TREE_COLORS\nSEPARATOR COMMA\nLEGEND_TITLE,Family\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_fam_legend.txt")
write.table(itol_fam2,"analysis/09_rhizobiales_phylogeny/iqtree/itol_fam_legend.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### nifH

nifh_df<-lables_both %>% mutate(nifh=ifelse(results==1,1,-1)) %>%
  select(new_label,nifh)
nifh_df$new_label<-str_replace_all(nifh_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,NifH\nCOLOR,#7ca118\nFIELD_SHAPES,2\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_nifh.txt")
write.table(nifh_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_nifh.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### source
source_df<-lables_both %>% select(new_label,source) %>% mutate(type="label",color=ifelse(source=="own","#ff0000","#000000")) %>% select(-source)
source_df$new_label<-str_replace_all(source_df$new_label," ","_")

cat("TREE_COLORS\nSEPARATOR COMMA\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_source.txt")
write.table(source_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_source.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### methane monooxygenase pmoC
pmoc_df<-lables_both %>% mutate(pmoc2=ifelse(pmoc==1,1,-1)) %>%
  dplyr::select(new_label,pmoc2)
pmoc_df$new_label<-str_replace_all(pmoc_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,particulate methane monooxygenase pmoC\nCOLOR,#eb4034\nFIELD_SHAPES,1\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_pmoc.txt")
write.table(pmoc_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_pmoc.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### methane monooxygenase mmoX
mmox_df<-lables_both %>% mutate(mmox2=ifelse(mmox==1,1,-1)) %>%
  dplyr::select(new_label,mmox2)
mmox_df$new_label<-str_replace_all(mmox_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,soluble methane monooxygenase mmoX\nCOLOR,#8e0e6a\nFIELD_SHAPES,1\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_mmox.txt")
write.table(mmox_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_mmox.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### methanol dehydrogenase xoxF
xoxf_df<-lables_both %>% mutate(xoxf2=ifelse(xoxf==1,1,-1)) %>%
  dplyr::select(new_label,xoxf2)
xoxf_df$new_label<-str_replace_all(xoxf_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,methanol dehydrogenase xoxF\nCOLOR,#ffdbac\nFIELD_SHAPES,3\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_xoxf.txt")
write.table(xoxf_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_xoxf.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)

### methanol dehydrogenase mxaF
mxaf_df<-lables_both %>% mutate(mxaf2=ifelse(mxaf==1,1,-1)) %>%
  dplyr::select(new_label,mxaf2)
mxaf_df$new_label<-str_replace_all(mxaf_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,methanol dehydrogenase mxaF\nCOLOR,#783f04\nFIELD_SHAPES,3\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_mxaf.txt")
write.table(mxaf_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_mxaf.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)

### gamma-glutamylmethylamide synthetase gmaS
gmas_df<-lables_both %>% mutate(gmas2=ifelse(gmas==1,1,-1)) %>%
  dplyr::select(new_label,gmas2)
gmas_df$new_label<-str_replace_all(gmas_df$new_label," ","_")

cat("DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,gamma-glutamylmethylamide synthetase gmaS\nCOLOR,#3d85c6\nFIELD_SHAPES,5\nFIELD_LABELS,f1\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_gmas.txt")
write.table(gmas_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_gmas.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


### number of metagenome occurrences, for our MAGs

mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mags_role$bac_phylum <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="p")


occurrences<-mags_role %>% filter(breadth>=50,Genome %in% gtdbtk$Genome) %>%
  group_by(Genome) %>% summarize(occurrences=n()) 

occurrences_df<-lables_both %>% left_join(occurrences) %>% select(new_label,occurrences)
occurrences_df$new_label<-str_replace_all(occurrences_df$new_label," ","_")


cat("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Occurrences\nCOLOR,#ff0000\nDATA\n",file="analysis/09_rhizobiales_phylogeny/iqtree/itol_occurrences.txt")
write.table(occurrences_df,"analysis/09_rhizobiales_phylogeny/iqtree/itol_occurrences.txt",append=TRUE,sep=",",quote = F, row.names = F, col.names=F)


