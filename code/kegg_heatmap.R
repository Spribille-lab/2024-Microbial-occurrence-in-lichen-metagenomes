#heatmap of proteins of interest

#setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(patchwork)
library(ape)
library(treeio)
require("phytools")
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(RColorBrewer)
options(repr.plot.width=20, repr.plot.height=40)
theme_set(theme_minimal(base_size = 23))

# 1. load data

##mag info
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())
annotated_mags_arranged<-annotated_mags %>% 
  group_by(bac_family) %>%
  arrange(bac_genus2, .by_group=T)


## kegg
kegg_of_interest<-read.delim("analysis/07_annotate_MAGs/kegg_of_interest.txt")
kegg_combined<-read.delim("analysis/07_annotate_MAGs/summarized_outputs/kegg_combined.txt")

## FeGenie annotation
fegenie<-read.csv("analysis/07_annotate_MAGs/fegenie/arkadiy/FeGenie-heatmap-data.csv")
fegenie$fegenie_class<-row.names(fegenie)

colnames(fegenie)[1]<-"fegenie_class"
fegenie$fegenie_class<-row.names(fegenie)

fegenie<-fegenie %>% pivot_longer(-fegenie_class,values_to="n",names_to="fasta")
fegenie$Genome<-fegenie$fasta %>% str_replace(".faa","")
fegenie <- fegenie %>% select(-fasta) %>% 
  pivot_wider(names_from = fegenie_class, values_from = n, values_fill=0 )
  

## balstp annotation
blastp<-read.delim("analysis/07_annotate_MAGs/summarized_outputs/manual_blast_search_results.txt")
blastp<-blastp[,1:5]
blastp<-blastp %>% filter(Genome %in% annotated_mags$Genome)

# 2. Combine data
##make a list of KOs, excluding associated with  iron, since for that, we have fegenie output
#kegg_list<-kegg_of_interest[kegg_of_interest$module_informal!="Iron transporter",1]
kegg_list<-kegg_of_interest[,1]


#make a table for all genes used here for annotation
fegenie_info<-data.frame(curr_name=colnames(fegenie[-1]),final_name=colnames(fegenie[-1]),type="iron, output from FeGenie")
kegg_info<-kegg_of_interest %>% #filter(module_informal!="Iron transporter") %>%
  mutate(curr_name=KO,final_name=paste(KO,gene,module_informal,sep=" ")) %>%
  select(curr_name,final_name,type)
blast_info<- blastp %>% mutate(curr_name=name,final_name=name)  %>%
  select(curr_name,final_name,type) %>% distinct()

gene_info<-rbind(kegg_info,fegenie_info,blast_info)

colnames(kegg_combined)<-c("locus","name","Genome")

df<-kegg_combined %>% filter(name %in% kegg_info$curr_name) %>%
  rbind(blastp %>% select(locus,name,Genome)) %>%
   group_by(Genome,name) %>%summarize(n=n()) %>% 
  pivot_wider(names_from=name,values_from=n,values_fill=0) %>%
  left_join(fegenie) 



# 3. Visualize
##prepare matrix
M <- t(as.matrix(df[,2:ncol(df)]))
colnames(M) = df$Genome
new_rownames<-data.frame(rownames(M)) %>% left_join(gene_info,by=c("rownames.M."="curr_name"))
rownames(M)<-new_rownames$final_name
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
type<-data.frame(rownames(M)) %>% left_join(gene_info,by=c("rownames.M."="final_name"))
type<-type$type
type_ordered<- gene_info$type %>% unique()

type_color<-c( "carbon metabolism"  = brewer.pal(8,"Pastel2")[1],
                      "photosynthesis"  = brewer.pal(8,"Pastel2")[2],
                      "nitrogen metabolism"  = brewer.pal(8,"Pastel2")[3],
                      "cofactors"  = brewer.pal(8,"Pastel2")[4],
                      "sulfur metabolism"  = brewer.pal(8,"Pastel2")[5],
                      "transporters"  = brewer.pal(8,"Pastel2")[6],
                      "iron, output from FeGenie"  = brewer.pal(8,"Pastel2")[7])
                      
ra = rowAnnotation(df=data.frame(type = type),
                  col = list(type = type_color),
                   annotation_legend_param = list(type=list(labels=type_ordered,at=type_ordered,title = "function")))


### set color scheme
tmp<-df %>% pivot_longer(-Genome,names_to="family",values_to="count")
col_fun = colorRamp2(c(0,1, max(tmp$count)), c("white","orange", "red"))


hm<-Heatmap(M, show_column_names = F, name = " ",
            row_order = gene_info$final_name[gene_info$final_name %in% rownames(M)],
            column_order=annotated_mags_arranged$Genome,
            top_annotation = ta,
            right_annotation = ra,
            col = col_fun,
            row_names_side = "left",
            heatmap_legend_param = list(
              title = "# of Genes", at = c(0, 1, 15, 30, 45), 
              labels = c(0, 1, 15, 30, 45))
)
hm
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/kegg_heatmap.pdf",width=20,height=30)
draw(hm,annotation_legend_side = "top")
dev.off()

## save human-readable table
df %>% pivot_longer(-Genome,names_to="family",values_to="count") %>%
  pivot_wider(names_from = Genome,values_from=count) %>%
  left_join(gene_info,by=c("family"="curr_name"))

#human-readable table
hr_table<-cazy_summarized %>% pivot_wider(names_from=family,values_from=n,values_fill = 0)  %>% 
  left_join(annotated_mags) 

#order columns and rows
hr_table<-hr_table %>%  select(sort(current_vars())) %>% 
  relocate(Genome,bac_genus2,bac_family,bac_order) %>% group_by(bac_family) %>%
  arrange(bac_genus2, .by_group=T)
write.table( hr_table, "analysis/07_annotate_MAGs/summarized_outputs/cazymes_summarized.txt", sep='\t',quote = F, row.names = F, col.names = T)

