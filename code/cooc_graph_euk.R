# coocurrence euks + cyano
#loading required packages
# All new packages should be added here
#setwd("~/Documents/coverage")
suppressPackageStartupMessages(require(tidyverse))
library(igraph)
library(qgraph)

# default plotting settings
options(repr.plot.width=10, repr.plot.height=5)
theme_set(theme_minimal(base_size = 23))


# define colors for plotting
node_colors = c("Lichen Fungal Symbiont" = "#FFC61E",
  "Green Algal Photobiont" = "#00CD6C", 
                "Cyanobacteria"="#09b1db",
                "Cyphobasidium"="#c90076",
                "Tremella"="#e091bf"
)

#load data
# here we load all data
# Manually curated table with functional assignments of the genomes (mycobiont/photobiont/etc) (produced by code/assigne_putative_mag_roles.R)
mags_role<-read.delim("analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.txt")
mags_role$label<-mags_role$confirmed_role
mags_role$label[mags_role$label=="photobiont_chloro"]<-"Green Algal Photobiont"
mags_role$label[mags_role$label=="photobiont_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$label=="cephalodia_cyano"]<-"Cyanobacteria"
mags_role$label[mags_role$label=="mycobiont"]<-"Lichen Fungal Symbiont"
mags_role$label[mags_role$label=="cypho"]<-"Cyphobasidium"
mags_role$label[mags_role$label=="tremella"]<-"Tremella"

#remove "fungi_other", as some of them are mycobionts in other metagenomes, and including them will be a misrepresentation
#remove metagenomes which have "mycobiont_missassigned" label, because they likely are not the lichen they are labelled as
#remove metagenomes where mycobiont mag wasn't recovered
mags_role_filtered<-mags_role %>% filter(!(label %in% c("fungi_other","algae_other","bacteria_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$label=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Lichen Fungal Symbiont"])

# clades to plot
clades_to_plot = c("Green Algal Photobiont", "Lichen Fungal Symbiont","Cyphobasidium","Tremella", "Cyanobacteria")



# Create table that collects which MAG is co-occurring with what other MAG
dat = mags_role_filtered %>% as_tibble() %>% select(Genome,metagenome)


# Some R Magic to create a who with who table
V <- crossprod(table(dat[c(2,1)]))
V_orig = V # save this for later, when we need the full object again
# remove duplicates and diagonal
V[upper.tri(V) ]<- 0
diag(V) <- 0


# for a graph we create edges 
edges = V %>% as_tibble() %>%
  mutate(from = rownames(V)) %>%
  gather(to, width, -from) %>%
  filter(width > 0) 

# lets define interesting groups for the graph annotation
species_labels = tibble(Genome = unique(c(edges$from, edges$to))) %>% 
  left_join(bind_rows(mags_role_filtered %>% select(Genome, label))) %>% unique()
# and nodes
vertices = tibble(name = unique(c(edges$from, edges$to))) %>% 
  left_join(species_labels, by = c("name" = "Genome")) %>%
  dplyr::rename(group = label) %>%
  mutate(color = node_colors[match(group, names(node_colors))])


# reduce edges and vertices to the predifined groups
vertices = vertices %>% filter(group %in% clades_to_plot)
edges = edges %>% filter(to %in% vertices$name & from %in%  vertices$name )



g = graph_from_data_frame(edges,vertices = vertices, directed = F)
# remove isoated verices
g = delete.vertices(g, which(degree(g)==0))

options(repr.plot.width=15, repr.plot.height=10)


# lets plot the graph nicely using qgraph
e <- get.edgelist(g, names=F)
# Make a layout that looks good
l = qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                      area=30*(vcount(g)^2),repulse.rad=(vcount(g)^3.6))


svg("results/figures/coocc_graph_euk.svg")
plot<-plot(g,layout=l,vertex.size=4,vertex.label=NA, weight=E(g)$weight)
op <- par(cex = 1.5)
legend(-1.6,1.5,, legend=names(node_colors),box.lty=0,bg=NA,
       fill=node_colors, cex=1)
text(0,-1.2,"Eukaryotic symbionts and Cyanobacteria")
dev.off()

##what are the two co-occurring algae?
edges[which.max(edges$width),]
