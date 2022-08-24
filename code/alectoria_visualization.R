#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(tidyverse)
library(plotly)
source('code/utils.R')

#paths to the tables
blast_path_X12<-"analysis/01_subsampling/X12/reports/blast_report.txt"
bbduk_path_X12<-"analysis/01_subsampling/X12/reports/bbduk_report.txt"
quast_path_X12<-"analysis/01_subsampling/X12/reports/quast_report.txt"

blast_path_GT57<-"analysis/01_subsampling/GT57/reports/blast_report.txt"
bbduk_path_GT57<-"analysis/01_subsampling/GT57/reports/bbduk_report.txt"
quast_path_GT57<-"analysis/01_subsampling/GT57/reports/quast_report.txt"

#get data frame
df_X12<-compile_data(blast_path_X12,bbduk_path_X12,quast_path_X12)
df_X12$metagenome<-"X12"
df_GT57<-compile_data(blast_path_GT57,bbduk_path_GT57,quast_path_GT57)
df_GT57$metagenome<-"GT57"

df<-rbind(df_X12,df_GT57)


#summarize
summary_wide<-df %>% group_by(yeast,bp,metagenome) %>% summarise(read_detection=mean(read_abs),blast_detection=mean( blast_abs),
        cypho_mean_coverage=mean(cypho_genome),tremella_mean_coverage=mean(tremella_genome))


#summary_detection
summary_wide<-df %>% group_by(yeast,bp,metagenome) %>% summarise(read_detection=mean(read_abs),blast_detection=mean( blast_abs))
summary_long<-summary_wide %>% pivot_longer(-c(yeast,bp,metagenome),names_to="variable",values_to="value")
summary_long$bp_char<-as.factor(summary_long$bp)

#heatmap
#function to make labels as percent
pct = function(x, digits=1) {
  sprintf(paste0("%1.", digits, "f%%"), x*100)
}

#fix facet labels
yeast.label<- c("Cyphobasidium, rDNA", "Tremella, rDNA")
names(yeast.label) <- c("cypho", "tremella")

detection.label <- c("Detection\nin assemblies","Detection\nin reads")
names(detection.label) <- c("blast_detection", "read_detection")

summary_long$metagenome<-factor(summary_long$metagenome, levels = c( "GT57","X12"), 
                      labels = c( "Cortex Slurry\nmetagenome","Bulk-lichen\nmetagenome"))

# Heatmap
color_vector<-c(col_GT57,col_X12,col_GT57,col_X12)
heatmap<-ggplot(summary_long,aes(x=bp,y=metagenome,fill=yeast))+
  geom_tile(fill="white")+geom_tile(aes(alpha=value),color="grey")+
  geom_text(aes(label=pct(value, 0)),size=3)+
facet_grid(variable~yeast,labeller = labeller(variable = detection.label, yeast = yeast.label))+
    scale_alpha_continuous(range=c(0.1,1), guide = FALSE)+
  scale_fill_manual(values=c(col_cypho,col_tremella),guide=F)+
  scale_x_log10(breaks = unique(summary_long$bp),labels=c("25\nMbp","50\nMbp","100\nMbp","250\nMbp","500\nMbp","1\nGbp","2.5\nGbp","5\nGbp","10\nGbp"),limits=c(16500000,15000000000))+
theme_minimal()+theme(strip.text.x = element_text(size = 12),
                      strip.text.y = element_text(size = 12,angle = 0),
                      plot.title = element_text(size=15),
                      axis.text.y = element_text(hjust = 1, colour = color_vector,size=12),
                      axis.title.y = element_blank(),
                      axis.title.x = element_blank(),
                      plot.margin=unit(c(0.3,0.3,0.5,0.3),"cm"))+
  ggtitle("A")

#summary coverage
summary_cov<-df %>% group_by(bp,metagenome) %>% summarise(cypho_mean_coverage=mean(cypho_genome),tremella_mean_coverage=mean(tremella_genome))
summary_cov_long<-summary_cov %>% pivot_longer(c(cypho_mean_coverage,tremella_mean_coverage),names_to="yeast",values_to="coverage")
summary_cov_long$yeast<-factor(summary_cov_long$yeast, levels = c( "cypho_mean_coverage","tremella_mean_coverage"), 
                                labels = c( "Cyphobasidium","Tremella"))

#plot coverage
coverage<-ggplot(summary_cov_long,aes(x=bp,y=coverage,color=metagenome))+
  geom_line()+geom_point()+
  scale_x_log10(breaks = unique(summary_cov$bp),labels=c("25\nMbp","50\nMbp","100\nMbp","250\nMbp","500\nMbp","1\nGbp","2.5\nGbp","5\nGbp","10\nGbp"),limits=c(16500000,15000000000))+
  scale_color_manual(values=c(col_GT57,col_X12))+
  scale_y_continuous(limits = c(0,40))+
  theme_minimal()+theme(legend.position = "none",
                        plot.title = element_text(size=15),
                        strip.text.x = element_text(size=12,face = "italic"),
                        plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +
  ylab("Genome coverage")+xlab("Sequencing depth, bp")+
  facet_wrap(yeast~.)+ggtitle("B")
  

#add annotations
ann_text1 <- data.frame(yeast="Cyphobasidium",bp=3000000000,coverage=37,metagenome="GT57",lab = "Cortex Slurry\nmetagenome")
ann_text2 <- data.frame(yeast="Cyphobasidium",bp=5000000000,coverage=7,metagenome="X12",lab = "Bulk-lichen\nmetagenome")
#ann_text3 <- data.frame(yeast="Tremella",bp=4000000000,coverage=13,metagenome="GT57",lab = "Cortex Slurry\nmetagenome")
#ann_text4 <- data.frame(yeast="Tremella",bp=5000000000,coverage=4,metagenome="X12",lab = "Bulk-lichen\nmetagenome")
#ann_text<-rbind(ann_text1,ann_text2,ann_text3,ann_text4)
ann_text<-rbind(ann_text1,ann_text2)
coverage<-coverage+geom_text(data = ann_text,label = ann_text$lab)

#combine plots
fig2<-heatmap /  coverage

ggsave("results/figures/alectoria_vis.svg",width=10,height=6)
ggsave("results/figures/alectoria_vis.png",width=10,height=6)
