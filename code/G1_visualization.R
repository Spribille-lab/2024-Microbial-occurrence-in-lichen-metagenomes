#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(gridExtra)
library(patchwork)
source('code/utils.R')

#get data set
blast_path<-"analysis/02_G1_visualization/blast_report_20210405.txt"
bbduk_path<-"analysis/02_G1_visualization/bbduk_report_26022021.txt"
quast_path<-"analysis/02_G1_visualization/quast_report_20210405.txt"

df<-compile_data(blast_path,bbduk_path,quast_path)


df2 <-df %>% group_by(yeast,bp) %>% summarize(mean_lec_cov=mean(lecanoro_genome),detection_blast=mean(blast_abs),
                                                    detection_reads=mean(read_abs),mean_cypho_cov=mean(cypho_genome),mean_alga_cov=mean(alga_genome))
df_fig1<-df2 %>% pivot_longer(c(detection_reads,detection_blast),names_to='method',values_to='detection')




#plot

#detection its
  df_fig1$yeast<-factor(df_fig1$yeast, levels = c("cypho", "tremella"), 
                        labels = c("Cyphobasidium", "Tremella"))
  
  df_fig1_split<-split(df_fig1,df_fig1$yeast)
 detection_cypho<- ggplot(data=df_fig1_split$Cyphobasidium,aes(x=bp,y=detection))+
    geom_line(aes(alpha=method),show.legend = F,color=col_cypho)+
    geom_area(aes(alpha=method), position = "identity",fill=col_cypho)+
    scale_x_log10(breaks = unique(df_fig1$bp),labels=c("25 Mbp","50 Mbp","100 Mbp","250 Mbp","500 Mbp","1 Gbp","2.5 Gbp","5 Gbp","10 Gbp"))+
    scale_alpha_manual(labels=c("In assemblies","In reads"),values=c(1,0.6))+
   scale_y_continuous(labels = scales::percent)+
    theme_minimal() +ggtitle("Detection by rDNA, Cyphobasidium")+theme(plot.title = element_text(hjust = 0.5))+
   ylab("Detection rate")+xlab("")
           
 detection_tremella<- ggplot(data=df_fig1_split$Tremella,aes(x=bp,y=detection))+
   geom_line(aes(alpha=method),show.legend = F,color=col_tremella)+
   geom_area(aes(alpha=method), position = "identity",fill=col_tremella)+
   scale_x_log10(breaks = unique(df_fig1$bp),labels=c("25 Mbp","50 Mbp","100 Mbp","250 Mbp","500 Mbp","1 Gbp","2.5 Gbp","5 Gbp","10 Gbp"))+
   scale_alpha_manual(labels=c("In assemblies","In reads"),values=c(1,0.6))+
   scale_y_continuous(labels = scales::percent)+
   ggtitle("Detection by rDNA, Tremella")+
   theme_minimal() + theme(plot.title = element_text(hjust = 0.5))+
   ylab("Detection rate")+xlab("Sequencing depth, bp")
 
#genome coverage
 df_fig2<-df2 %>%filter(yeast=="cypho") %>% pivot_longer(c(mean_alga_cov,mean_cypho_cov,mean_lec_cov),names_to="genome",values_to="coverage")
 df_fig2$genome<-factor(df_fig2$genome, levels = c("mean_lec_cov","mean_alga_cov","mean_cypho_cov"), 
                       labels = c("Lecanoromycete","Alga","Cyphobasidium"))
 
 coverage<-ggplot(data=df_fig2)+geom_line(aes(x=bp,y=coverage,color=genome))+
   geom_text(data=subset(df_fig2, coverage>0),aes(label = ifelse(coverage < 0, NA, coverage), x=bp, y = coverage,color=genome), position = "stack",show.legend = F,vjust=-.6)+
  scale_x_log10(breaks = unique(df_fig1$bp),labels=c("25 Mbp","50 Mbp","100 Mbp","250 Mbp","500 Mbp","1 Gbp","2.5 Gbp","5 Gbp","10 Gbp"))+
  scale_color_manual(values=c(col_lecanoro,col_alga,col_cypho)) +
   ggtitle("Average Genome Coverage")+
   theme_minimal() + theme(plot.title = element_text(hjust = 0.5))+
   ylab("Average Coverage")+xlab("Sequencing depth, bp")
   
#assemble plots           
 detection<-detection_cypho /  detection_tremella
 detection / coverage + plot_layout(heights = c(1, 1,3))
 
ggsave("results/figures/G1_vis.svg",height=10)
ggsave("results/figures/G1_vis.png",height=10)


