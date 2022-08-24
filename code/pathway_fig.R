#heatmap of proteins of interest

setwd("~/Documents/coverage")
source("code/utils.R")
library(tidyverse)
library(RColorBrewer)
library(svglite)
options(repr.plot.width=20, repr.plot.height=40)
theme_set(theme_minimal(base_size = 23))

# 1. load mag info


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
#remove non-complete mags added out of curiosity
annotated_mags_arranged<- annotated_mags_arranged %>% filter(Genome!="public_SRR14722130_metawrap_bin.2" & Genome!="private_T1889_metawrap_bin.7")

#2. kegg annotations
###load data
kegg_of_interest<-read.delim("analysis/07_annotate_MAGs/kegg_of_interest_fig.txt")
kegg_combined<-read.delim("analysis/07_annotate_MAGs/summarized_outputs/kegg_combined.txt")
colnames(kegg_combined)<-c("locustag","KO","Genome")

###prepare dataset
kegg_list<-kegg_of_interest[,1]

kegg_df <- kegg_combined %>% filter(KO %in% kegg_list) %>% group_by(Genome,KO) %>%
  summarize(sum=n()) %>% mutate(presence=ifelse(sum>0,1,0)) %>% dplyr::select(-sum) %>%
  pivot_wider(names_from = KO,values_from = presence,values_fill = 0) %>%
#new variable: calvin cycle
    mutate(Calvin_cycle = ifelse(K00855==1 & K01601 == 1 & K01602 ==1 & K00927 == 1 &
             K00134==1 & (K01623==1 | K01624==1) & K00615 ==1 & K03841 ==1 &
               (K01807==1 | K01808==1),1,0)) %>%
  #new variable: anoxygenic photosystem II
  mutate(photosystem = ifelse(K08928==1 & K08929==1 & K13991==1 & K13992==1 &
                                K08926==1 & K08927==1,1,0)) %>%
  #new variable: bacteriochlorophyll
  mutate(bacteriochlorophyll = ifelse(K04035==1 & K04037==1 & K04038==1 & K04039==1 &
         K11333==1 & K11334==1 & K11335==1 & K11336==1 & K11337==1 & K04040==1 & K10960==1,1,
        ifelse(K04035==1 & K04037==1 & K04038==1 & K04039==1 &
                  K11334==1 & K11335==1 & K11336==1 & K11337==1 & K04040==1 & K10960==1,0.9,0))) %>%
  #new variable: carotenoid
  mutate(carotenoids = ifelse(K02291==1 & K10027==1 & K09844==1 & 
         K09844==1 & K09845==1 & K09846==1,1,0)) %>%
  #new variable: nitrogen fixation
  mutate(nitrogen_fixation = ifelse(K02588==1,1,0)) %>%
  #new variable: Cobalamin
  mutate(cobalamin = ifelse((K00768==1 & K02226==1) | 
         (K02232==1 & K02231==1 & (K02227==1 | K02225==1) & (K00798==1 | K19221==1) ),
     1,ifelse(K02232==1 & K02231==1  & (K00798==1 | K19221==1),0.9,0)
     )) %>%
  #new variable: biotin
  mutate(biotin = ifelse(K00652+K00833+K01935+K01012==4,1,
                         ifelse(K00652+K00833+K01935+K01012==3,0.75,0))) %>%
  #new variable: Riboflavin
  mutate(riboflavin = ifelse((K00794==1 & K00793==1 & K11753==1) |
         ((K01497==1 | K14652==1) & K11752==1 & K21064==1),1,0)) %>%
  #new variable: thiamin
  mutate(thiamine = ifelse(K00878==1 & K00941==1 & K00788==1,1,0)) %>%
  #new variable: sorbitol/mannitol transporter
  mutate(sorbitol_mannitol_transporter = ifelse(K10227==1 & 
         K10228==1 & K10229==1 & K10111==1,1,0)) %>%
  #new variable: urea transporter
  mutate(urea_transporter = ifelse(K11959==1 &
                                   K11960==1 &
                                   K11961==1 &
                                   K11962==1 &
                                   K11963==1,1,0)) %>%
  #new variable: Erythritol transporter
  mutate(erythritol_transporter = ifelse(K17202==1 &
                                         K17203==1 &
                                         K17204==1,1,0)) %>%
  #new variable: xylitol transporter
  mutate(xylitol_transporter = ifelse(K17205==1 &
                                      K17206==1 &
                                      K17207==1,1,0)) %>%
  #new variable: Inositol  transporter
  mutate(inositol_transporter = ifelse(K17208==1 &
                                       K17209==1 &
                                       K17210==1,1,0)) %>%
  
  #new variable: glycerol  transporter
  mutate(glycerol_transporter = ifelse(K17321==1 &
                                       K17322==1 &
                                       K17323==1 &
                                       K17324==1 &
                                       K17325==1,1,0)) %>%
  #new variable: urease
  mutate(urease = ifelse((K14048==1 | (K01430==1 & K01429==1)) & K01428==1,1,0)) %>%
  #new variable: Fucose transporter
  mutate(fucose_transporter = ifelse(K02429==1,1,0)) %>%
  #new variable: Glycerol aquaporin transporter
  mutate(glycerol_aquaporin_transporter = ifelse(K02440==1,1,0)) %>%
  #new variable: Glycerol/sorbitol transporter
  mutate(glycerol_sorbitol_transporter = ifelse(K02781==1 &
                                                K02782==1 &
                                                K02783==1,1,0)) %>%
  #new variable: ammonium transporter
  mutate(ammonium_transporter = ifelse(K03320==1,1,0)) %>%
  #new variable: ribose transporter
  mutate(ribose_transporter = ifelse(K10439==1 &
                                     K10440==1 &
                                     K10441==1,1,0)) %>%
  #new variable: xylose transporter
  mutate(xylose_transporter = ifelse(K10543==1 &
                                     K10544==1 &
                                     K10545==1,1,0)) %>%
  #new variable: multiple sugar transporter
  mutate(multiple_sugar_transporter = ifelse(K10546==1 &
                                             K10547==1 &
                                             K10548==1,1,0)) %>%
  #new variable: fructose transporter
  mutate(fructose_transporter = ifelse(K10552==1 &
                                       K10553==1 &
                                       K10554==1,1,0)) %>%
  #new variable: arabinose transporter
  mutate(arabinose_transporter = ifelse(K10537==1 &
                                        K10538==1 &
                                        K10539==1,1,0)) %>%
  #new variable: branched-chain amino acid transporter
  mutate(branched_transporter = ifelse(K01999==1 &
                                       K01997==1 &
                                       K01998==1 &
                                       K01995==1 &
                                       K01996==1,1,0)) %>%
  #new variable: L-amino acid transporter
  mutate(l_amino_transporter = ifelse(K09969==1 &
                                      K09970==1 &
                                      K09971==1 &
                                      K09972==1,1,0)) %>%
  #new variable: glutamate transporter
  mutate(glutamate_transporter = ifelse(K10001==1 &
                                        K10002==1 &
                                        K10003==1 &
                                        K10004==1,1,0)) %>%
  #new variable: capsular transporter
  mutate(capsular_transporter = ifelse(K10107==1 &
                                       K09688==1 &
                                       K09689==1,1,0)) %>%
  #new variable: methylotrophy
  mutate(methanol_dehydrogenase = ifelse(K23995==1,1,0)) %>%
  #remove KO columns
  dplyr::select(-contains("K"))# %>% left_join( annotated_mags_arranged)
  
  

# 4.  viz



df_long<-kegg_df %>% pivot_longer(-Genome, names_to = "pathway",values_to = "presence")

df_long<- annotated_mags_arranged %>% 
  inner_join(df_long) %>% mutate(presence_factor = ifelse(
  presence==1,"full",ifelse(presence>0,"partial","missing"))) %>%
  mutate(presence_size = ifelse(
    presence>0,1,NA)) %>%
  mutate(bac_family_label = ifelse(bac_family=="Sphingomonadaceae","Sphingo\nmonadaceae",bac_family)) 
  
#order rows and columns
df_long$Genome<-factor(df_long$Genome,level = annotated_mags_arranged$Genome)

pathway_type<-data.frame("pathway"=c("capsular_transporter",
                                     "branched_transporter" ,         
                                     "l_amino_transporter",           
                                     "glutamate_transporter" ,
                                     "ammonium_transporter",
                                     "urea_transporter",
                                     "ribose_transporter" ,           
                                     "xylose_transporter" ,           
                                     "multiple_sugar_transporter"  ,  
                                     "fructose_transporter" ,         
                                     "arabinose_transporter" ,
                                     "fucose_transporter",
                                     "erythritol_transporter",        
                                     "xylitol_transporter",           
                                     "inositol_transporter",          
                                     "glycerol_transporter",
                                     "glycerol_aquaporin_transporter",
                                     "glycerol_sorbitol_transporter",                      
                                     "sorbitol_mannitol_transporter" ,              
                                     "thiamine",
                                     "riboflavin"  ,   
                                     "cobalamin",
                                     "biotin" ,
                                     "urease",
                                     "nitrogen_fixation",
                                     "methanol_dehydrogenase",
                                     "Calvin_cycle",                  
                                     "carotenoids",
                                     "bacteriochlorophyll",
                                     "photosystem"
),"functions"=c(rep("Other transporters",6),rep("Carbohydrate transporters",13),
                rep("Cofactors",4),rep("C and N metabolism",4), rep("Photo-\nsynthesis",3)
))

df_long<-df_long%>% left_join(pathway_type)
df_long$pathway<-factor(df_long$pathway,level = c("capsular_transporter",
                                                  "branched_transporter" ,         
                                                          "l_amino_transporter",           
                                                          "glutamate_transporter" ,
                                                          "ammonium_transporter",
                                                          "urea_transporter",
                                                          "ribose_transporter" ,           
                                                          "xylose_transporter" ,           
                                                          "multiple_sugar_transporter"  ,  
                                                          "fructose_transporter" ,         
                                                          "arabinose_transporter" ,
                                                          "fucose_transporter",
                                                          "erythritol_transporter",        
                                                          "xylitol_transporter",           
                                                          "inositol_transporter",          
                                                          "glycerol_transporter",
                                                          "glycerol_aquaporin_transporter",
                                                          "glycerol_sorbitol_transporter",                      
                                                          "sorbitol_mannitol_transporter" ,              
                                                         "thiamine",
                                                         "riboflavin"  ,   
                                                         "cobalamin",
                                                         "biotin" ,
                                                         "iron_ion_transport","siderophore_synthesis",
                                                          "urease",
                                                          "nitrogen_fixation",
                                                         "methanol_dehydrogenase",
                                                          "Calvin_cycle",                  
                                                          "carotenoids",
                                                          "bacteriochlorophyll",
                                                          "photosystem"
                                                          ))

df_long$functions<-factor(df_long$functions,level = c("Photo-\nsynthesis","C and N metabolism","Cofactors","Carbohydrate transporters", "Other transporters" ))

ggplot(df_long,aes(x=Genome,y=pathway,size=presence_size,shape=presence_factor,color=bac_family))+
  geom_point()+
   scale_shape_manual(values=c(16,10),limits = c("full","partial") )+
  guides(size="none",color="none",shape = guide_legend(override.aes = list(size = 6) ))+theme_minimal()+
  facet_grid(functions~bac_family_label,scales = "free",space="free",labeller = label_wrap_gen(width=5))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size =12),
        strip.text.y = element_text(size =12,angle = 0),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  xlab("")+ylab("")+
  labs(shape="Pathway\nCompleteness")+
  scale_y_discrete(labels=c("capsular_transporter"="Capsular polysaccharide Transporter kpsMTE",   
                            "branched_transporter" = "Branched-chain amino acid Transporter LivGFHKM" ,         
                            "l_amino_transporter"="L-amino acid Transporter AapJMPQ",           
                            "glutamate_transporter" ="Glutamate/Aspartate Transporter GltIJKL",
                            "ammonium_transporter"="Ammonium Transporter Amt",
                            "urea_transporter"="Urea Transporter urtABCDE",
                            "ribose_transporter" ="Ribose/D-xylose Transporter RbcABC",           
                            "xylose_transporter" ="D-xylose Transporter XylFGH",           
                            "multiple_sugar_transporter" ="Multiple sugar Transporter ChvE-GguAB" ,  
                            "fructose_transporter"="Fructose Transporter FrcABC" ,         
                            "arabinose_transporter"= "Arabinose Transporter AraFGH" ,
                            "fucose_transporter"="Fucose Transporter fucP",
                            "erythritol_transporter"="Erythritol Transporter EryEFG",        
                            "xylitol_transporter"="Xylitol Transporter XltABC",           
                            "inositol_transporter"="Inositol Transporter IatAP-IbpA",          
                            "glycerol_transporter" = "Glycerol Transporter GlpPQSTV",
                            "glycerol_aquaporin_transporter"="Glycerol Transporter GLPF",
                            "glycerol_sorbitol_transporter"="Glucitol/sorbitol Transporter SrlABE",                      
                            "sorbitol_mannitol_transporter"="Sorbitol/mannitol Transporter SmoEFGK" ,              
                            "thiamine"="Thiamine salvage pathway",
                            "riboflavin" ="Riboflavin biosynthesis" ,   
                            "cobalamin"="Cobalamin biosynthesis" ,
                            "biotin"="Biotin biosynthesis" ,
                            "urease"="Urease UreABC",
                            "nitrogen_fixation"="Nitrogenase NifH, Nitrogen fixation",
                            "methanol_dehydrogenase"="Methanol Dehydrogenase, Methylotrophy",
                            "Calvin_cycle"="Calvin cycle",                  
                            "carotenoids"="Carotenoid biosynthesis",
                            "bacteriochlorophyll"="Bacteriochlorophyll biosynthesis",
                            "photosystem"="Anoxygenic Phtosystem II"))

ggsave("results/figures/pathway_fig.svg",width=16,height=7,device="svg",bg="white")







