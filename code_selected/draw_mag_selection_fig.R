#figure to illustrate the process of MAG selection for annotation
library(tidyverse)
library(patchwork)
require(scales)
library(waffle)
library(extrafont)
library(hrbrthemes)

#1. load data

## number of occurrences for per mag
mag_occur<-read.delim("analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_mags_frequency.tsv")
colnames(mag_occur)[2]<-'mag_occ'

## number of occurrences for per genus
genus_occur<-read.delim("analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_genera_frequency.tsv")
genus_occur<- genus_occur %>% mutate(bac_genus2=ifelse(bac_genus2 =="Unknown gen. sp.",paste0(bac_order," fam. gen. sp."),bac_genus2))
colnames(genus_occur)[2]<-'genus_occ'
genus_occur$if_selected<-"not_selected" ##dummy variable to code whether the genus was selected
genus_occur$if_selected[1:13]<-"selected"

## MAG QC results
checkm<-read.delim("analysis/05_MAGs/tables/checkm_results.tab",header=F,col.names=c("Genome","compelteness","contamination","strain_heterog","taxonomy"))


#2. combine datasets together
df<-mag_occur  %>%
  left_join(genus_occur) %>%
  left_join(checkm)

#3. visualize: bars showing # of occurrences
#cols <- c("Beijerinckiaceae"="#D38E18",  "Acetobacteraceae"= "#CB1B06", "Sphingomonadaceae" = "#A37FBF", "Acidobacteriaceae" = "#2D4FB1", "Nostocaceae" = "#1DA6B5","UBA10450"="#8CB28D")

bar<-ggplot(genus_occur,aes(fill=if_selected,
                            y=reorder(bac_genus2, genus_occ),
                            x=genus_occ))+
  geom_col()+
  scale_fill_manual(name = "Selected",values=c("selected"="#990000","not_selected"="#bfbfbf"))+
  #scale_alpha_discrete(range=c(0.2, 1))+
  guides(fill = "none")+
  xlab("Number of occurrences per bacterial genus")+
  ylab("Genera")+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=7),
        )

##4. visualize: waffle chart showing selected MAGs
###4.1.loading fonts:followed https://rud.is/rpubs/building-pictograms.html
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#install_fa_fonts() copied cp  /Users/gulnaratagirdzhanova/Library/R/x86_64/4.1/library/waffle/fonts/fa-solid-900.ttf /Library/Fonts
font_import()
extrafont::loadfonts(quiet = TRUE)

###4.2. prepare new data frame
df <- df %>% mutate(completeness_binary=ifelse(compelteness>95,"complete","incomplete"))

df_summarized<-df %>% group_by(bac_genus2,completeness_binary) %>% summarize(n=n()) %>%
  filter(!is.na(bac_genus2)) %>% left_join(genus_occur) %>%
  mutate(mag_selected=if_else(if_selected=="selected"&completeness_binary=="complete","selected","not selected"))

#reorder 
library(forcats)
df_summarized$bac_genus2<-as.factor(df_summarized$bac_genus2)
df_summarized$bac_genus2<-reorder(df_summarized$bac_genus2, df_summarized$genus_occ)
df_summarized$bac_genus2<-fct_rev(df_summarized$bac_genus2)

### 4.3. draw
waffle<-ggplot(df_summarized %>% filter(bac_genus2 %in% genus_occur[1:18,1]),aes(label=completeness_binary,color=mag_selected,values=n))+
  waffle::geom_pictogram(n_rows = 2,flip=F,size=2, family = "Font Awesome 5 Free Solid")+
  scale_color_manual(values=c("selected"="#990000","not selected"="#bfbfbf"))+
  scale_label_pictogram(
    name = NULL,
    values = c("circle","circle-notch"),
    labels = c("complete","incomplete")
  ) +
  facet_wrap(bac_genus2~.,ncol = 1,strip.position = "left")+
  coord_equal() +
  guides(color = "none",label="none")+
  theme_enhance_waffle()+theme_minimal()+
  theme( panel.grid = element_blank(),
         axis.text= element_blank(),
         strip.text.y.left = element_text(angle = 0,hjust=1,size=7),
         legend.text = element_text(size=7))
  

bar+waffle + plot_layout(widths = c(1, 1))
ggsave("results/figures/mag_selection.svg",device="svg",width=180,height=120,unit="mm")






###exploratory
###1. vertical waffle bars
df_summarized<-df %>% group_by(bac_genus2,completeness_binary) %>% summarize(n=n()) %>%
  filter(!is.na(bac_genus2)) %>% left_join(genus_occur) %>%
  mutate(mag_selected=if_else(if_selected=="selected"&completeness_binary=="complete","selected","not selected")) %>%
  arrange(desc(completeness_binary))

#reorder 
df_summarized$bac_genus2<-as.factor(df_summarized$bac_genus2)
df_summarized$bac_genus2<-reorder(df_summarized$bac_genus2, df_summarized$genus_occ)
df_summarized$completeness_binary<-factor(df_summarized$completeness_binary,levels=c("incomplete","complete"))


ggplot(df_summarized %>% filter(bac_genus2 %in% genus_occur[1:25,1]),aes(label=completeness_binary,color=mag_selected,values=n))+
  waffle::geom_pictogram(n_rows = 3,flip=T,size=3)+
  scale_color_manual(values=c("selected"="#990000","not selected"="#bfbfbf"))+
  scale_label_pictogram(
    name = NULL,
    values = c("circle","circle-notch"),
    labels = c("complete","incomplete")
  ) +
  facet_wrap(~bac_genus2,nrow = 1, strip.position = "bottom")+
  theme_enhance_waffle()+
  theme( panel.grid = element_blank(),
         axis.text= element_blank())






