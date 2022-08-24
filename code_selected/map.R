#generate table for the map
setwd("~/Documents/coverage")
library(tidyverse)
library(stringr)
library(ggmap)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")


## 1. load data
mtg_info<-read.delim("results/tables/all_metagenome_reanalysis.txt")
location_raw<-read.delim("analysis/03_metagenome_reanalysis/locations_compiled_manually.txt")

### check that all metagenomes have coordinates
mtg_info %>% filter(!(Run %in% location_raw$run_accession))
location_raw<-location_raw %>% filter(run_accession %in%mtg_info$Run)

## 2. transform coord format: N/S/W/E format
###make a subset of coordinates in the wrong format
dataset1<-location_raw %>% filter(grepl("N",location)|grepl("S",location)) 

coord_tmp<-str_split(dataset1$location," ",simplify = TRUE) %>% data.frame()
coord_tmp<-cbind(coord_tmp,dataset1$run_accession)
colnames(coord_tmp)<-c("lat1","lat_type","long1","long_type","run_accession")

coord_tmp<-coord_tmp %>% mutate(lat=ifelse(lat_type=="N",as.numeric(lat1),-1*as.numeric(lat1))) %>%
  mutate(long=ifelse(long_type=="E",as.numeric(long1),-1*as.numeric(long1))) 

dataset1<-coord_tmp %>% select(run_accession,lat,long)

## 3. transform coord format: comma
dataset2 <- location_raw %>% filter(!(run_accession %in% dataset1$run_accession))

coord_tmp2<-str_split(dataset2$location,", ",simplify = TRUE) %>% data.frame()
coord_tmp2<-cbind(coord_tmp2,dataset2$run_accession)
colnames(coord_tmp2)<-c("lat","long","run_accession")

dataset2<-coord_tmp2

## 4. combine dataset
dataset<-rbind(dataset1,dataset2) %>% filter(long!="")
write.table(dataset,"analysis/03_metagenome_reanalysis/locations_final.txt",row.names = F,quote = F,sep="\t")


## 5. make a map
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf(fill= "antiquewhite",color="lightgrey",size=0.2)+
  coord_sf(expand = FALSE)+
  geom_point(data = dataset, aes(x = as.numeric(long), y = as.numeric(lat)), color = "brown", size = 2)+ 
  theme_minimal()+
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank(),
       axis.text.x = element_blank(),
       panel.background = element_rect(fill = "aliceblue"))

ggsave("results/figures/map.pdf",width=7.8,height=4.8)



li<-mags_role_filtered %>% filter(bac_genus2=="RH-AL1")
t<-dataset %>% filter(run_accession %in% li$metagenome)
ggplot(data = world) +
    geom_sf(fill= "antiquewhite",color="lightgrey",size=0.2)+
     coord_sf(expand = FALSE)+
     geom_point(data = dataset, aes(x = as.numeric(long), y = as.numeric(lat)), color = "brown", size = 2)+ 
     geom_point(data = t, aes(x = as.numeric(long), y = as.numeric(lat)), color = "red", size = 2)+ 
     theme_minimal()+
     theme(axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.x = element_blank(),
                     panel.background = element_rect(fill = "aliceblue"))
 View(li)

