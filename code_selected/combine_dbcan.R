#Creating combined kegg mapper files by genus

## 1. misc
library(tidyverse)
library(stringr)
library(seriation)
library(ComplexHeatmap)
library(DECIPHER)
library(circlize)
library(RColorBrewer)
source("code/utils.R")

## 2. read data
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
annotated_mags<-read.delim("analysis/07_annotate_MAGs/mag_list.txt",header=F,col.names = "Genome")

## 3. get taxonomy
mags_role$bac_genus <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g")
mags_role$bac_family <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="f")
mags_role$bac_order <- sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="o")
mags_role<-mags_role %>%mutate(bac_genus2=ifelse(bac_genus=="Unknown",paste(bac_family," gen. sp.",sep=""),bac_genus))

annotated_mags<-annotated_mags %>% left_join(mags_role %>% select(Genome, bac_genus2,bac_family,bac_order) %>% distinct())

## 4. make a function to process dbcan output

process_dbcan<-function(mag){

hmmer_filename<-paste0("analysis/07_annotate_MAGs/",mag,"_dbcan/hmmer.out")
hmmer<-read.delim(hmmer_filename)

hmmer_filtered<-hmmer %>% filter(E.Value<1e-20,Coverage>0.3)
hmmer_filtered$family<-str_replace(hmmer_filtered$HMM.Profile,".hmm","")

diamond_filename<-paste0("analysis/07_annotate_MAGs/",mag,"_dbcan/diamond.out")
diamond<-read.delim(diamond_filename)

diamond_filtered<-diamond %>% filter(E.Value<1e-20,X..Identical>30)
diamond_filtered<-diamond_filtered %>%  mutate(family = str_replace(CAZy.ID, "^(.+?)\\|","")) %>% 
  separate_rows(family,sep = "\\|", convert = FALSE)

combined<-diamond_filtered %>% select( Gene.ID, family,CAZy.ID) %>% inner_join(hmmer_filtered) %>%
  select(Gene.ID, family)
combined$Genome<-mag
return(combined)
}

## 5. apply the function to all mags
l<-lapply(annotated_mags$Genome,process_dbcan)
cazy_combined<-do.call(rbind,l)
cazy_combined$class<-substr(cazy_combined$family,1,2)
write.table(cazy_combined, "analysis/07_annotate_MAGs/summarized_outputs/cazymes_gene_assignments.txt", sep='\t',quote = F, row.names = F, col.names = T)

## 6. process the result
cazy_summarized<-cazy_combined %>% group_by(Genome,family) %>% summarize(n=n()) 

#human-readable table
hr_table<-cazy_summarized %>% pivot_wider(names_from=family,values_from=n,values_fill = 0)  %>% 
  left_join(annotated_mags) 

#order columns and rows
hr_table<-hr_table %>%  select(sort(current_vars())) %>% 
  relocate(Genome,bac_genus2,bac_family,bac_order) %>% group_by(bac_family) %>%
  arrange(bac_genus2, .by_group=T)
write.table( hr_table, "results/tables/cazymes_summarized.txt", sep='\t',quote = F, row.names = F, col.names = T)

# summarize by genus
family_genus<-cazy_combined %>% group_by(Genome,family) %>% summarize(n=n()) %>% left_join(annotated_mags) %>%
  group_by(bac_family,family) %>% summarize(mean=mean(n),n_mags=n())

# summarize by class and genus: calculate avergae # f cazymes in diff. classes, add median. # of cazymes total, and taxonomy
class<-cazy_combined %>% group_by(Genome,class) %>% summarize(n=n()) %>% left_join(annotated_mags)

class_genus <- class %>% group_by(bac_genus2,class) %>% 
  summarize(median=median(n)) 
total_cazy<-cazy_combined %>% group_by(Genome) %>%  summarize(n=n()) %>% left_join(annotated_mags) %>%
  group_by(bac_genus2) %>% summarize(total=median(n))

class_genus_summ <- class_genus %>% pivot_wider(values_from = median,names_from=class,values_fill=0) %>% left_join(total_cazy) %>% left_join(annotated_mags %>% select(-Genome) %>% distinct())
write.table( class_genus_summ, "results/tables/median_cazy_by_genus.txt", sep='\t',quote = F, row.names = F, col.names = T)

#in mags percantage of cazymes occupied by diff classes, median by family
total_per_mag<-class %>% group_by(Genome) %>% summarize(total=sum(n))
cazy_perc<-class %>% left_join(total_per_mag) %>% mutate(percent=n*100/total) %>% 
  group_by(bac_family,class) %>% summarize(avg_perc=median(percent)) %>%
  pivot_wider(values_from=avg_perc,values_fill=0,names_from=class)
write.table( cazy_perc, "analysis/07_annotate_MAGs/summarized_outputs/cazy_class_percentage_by_bac_family.txt", sep='\t',quote = F, row.names = F, col.names = T)


#exploratory
ggplot(class,aes(x=class,y=n)) + geom_point()+facet_wrap(~bac_genus2)
ggplot(class_genus %>% left_join(annotated_mags %>% select(-Genome) %>% distinct()),aes(x=class,y=mean)) + geom_point()+facet_wrap(~bac_family)

#saved: number of cazy from different classes in each mag grouped by family
ggplot(class,aes(x=class,y=n)) + geom_jitter(width=0.1)+facet_wrap(~bac_family)+
  xlab("")+ylab("Number of identified CAZymes")
ggsave("results/figures/cazy_bac_genus.png",device="png",height = 5,width=10,bg="white")


## 7. visualize: hm

###prepare matrix
cazy_hm<-cazy_summarized %>% pivot_wider(names_from=family,values_from=n,values_fill = 0)
M <- t(as.matrix(cazy_hm[,2:ncol(cazy_hm)]))
colnames(M) = cazy_hm$Genome

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


### make top annotation for the bacterial taxonomy
genus_df<-data.frame(colnames(M)) %>% left_join(annotated_mags %>%select(Genome, bac_genus2) %>%unique(),by=c("colnames.M."="Genome"))
genus<-genus_df$bac_genus2

family_df<-data.frame(colnames(M)) %>% left_join(annotated_mags %>%select(Genome, bac_family) %>%unique(),by=c("colnames.M."="Genome"))
family<-family_df$bac_family

genus_ordered<-hr_table %>%  ungroup() %>% select(bac_genus2) %>% distinct()
genus_ordered<-genus_ordered[,1]

ta = HeatmapAnnotation(df = data.frame(family=family,genus = genus), show_annotation_name=T,
         annotation_legend_param = list(genus=list(labels=genus_ordered,at=genus_ordered)),
                       col = list(genus = genus_node_color,family=family_node_color)
                       )

### set color scheme
col_fun = colorRamp2(c(0, max(cazy_summarized$n)), c("white", "red"))

###plot the heatmap
hm<-Heatmap(M, show_column_names = F, name = " ",
        row_order = colnames(hr_table)[-(1:4)],
        column_order=hr_table$Genome,
        top_annotation = ta,
        column_title = "              Fig. S8. Heatmap of CAZy families present in the bacterial MAGs. We selected highe-quality MAGs from the 13 most frequently occurring bacterial genera."
        #col = col_fun
        )
###export
pdf(file="analysis/07_annotate_MAGs/summarized_outputs/cazy_heatmap.pdf",width=15,height=40)
draw(hm)
dev.off()
