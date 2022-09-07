# coocurrence mycobiont + lichenihabitans
#loading required packages
# All new packages should be added here
#setwd("~/Documents/coverage")
suppressPackageStartupMessages(require(tidyverse))
library(igraph)
library(qgraph)
library(plotly)
source("code/utils.R")

# default plotting settings
options(repr.plot.width=10, repr.plot.height=5)
theme_set(theme_minimal(base_size = 23))


# define colors for plotting
node_colors = c("Main fungal symbiont" = "#FFC61E",
                "Lichenihabitans" = "#783f04"
)

#load data
# Manually curated table with functional assignments of the genomes (mycobiont/photobiont/etc) (produced by code/assigne_putative_mag_roles.R)
mags_role<-read.delim("results/tables/MAG_confirmed_roles_bwa.txt")
mags_role$label<-as.character(mags_role$confirmed_role)
mags_role$label[mags_role$label=="mycobiont"]<-"Main fungal symbiont"

#add label for licheninhabintans
mags_role$bac_genus = as.character(sapply(mags_role$bat_bacteria, gtdb_get_clade, clade="g"))
mags_role$label[mags_role$bac_genus=="Lichenihabitans"]<-"Lichenihabitans"


#remove "fungi_other", as some of them are mycobionts in other metagenomes, and including them will be a misrepresentation
#remove metagenomes which have "mycobiont_missassigned" label, because they likely are not the lichen they are labelled as
#remove metagenomes where mycobiont mag wasn't recovered
mags_role_filtered<-mags_role %>% filter(!(label %in% c("fungi_other","algae_other","bacteria_other"))) %>%
  filter(!(metagenome %in% mags_role$metagenome[mags_role$label=="mycobiont_missassigned"])) %>%
  filter(metagenome %in% mags_role$metagenome[mags_role$label=="Main fungal symbiont"])

# clades to plot
clades_to_plot = c("Lichenihabitans", "Main fungal symbiont")



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
  rename(group = label) %>%
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

pdf("results/figures/coocc_graph_lichenihab.pdf")
plot<-plot(g,layout=l,vertex.size=4,vertex.label=NA, weight=E(g)$weight)
op <- par(cex = 1.5)
legend(-1.6,1.5,, legend=names(node_colors),box.lty=0,bg=NA,
       fill=node_colors, cex=1)
text(0,-1.2,"Lichenihabitans (Beijerinckiaceae)")

dev.off()



### what are the three dominant licheninhabitans?
edges %>% arrange(desc(width))

edges %>% filter(from=="public_SRR14722130_metawrap_bin.2"&to=="private_T1889_metawrap_bin.7")

#how many nodes they connect to
edges %>% group_by(from) %>% summarize(n=n()) %>% arrange(desc(n))


genome_coverage = read_csv(paste0(basepath,"tables/read_mapping/bwa_coverage.csv")) 
genome_coverage2<-genome_coverage %>%pivot_longer(-Genome,values_to="coverage",names_to="metagenome") %>%
  pivot_wider(names_from=Genome,values_from = coverage) 

g1<-ggplot(genome_coverage2)+geom_point(aes(x=public_SRR14722130_metawrap_bin.2,y=private_T1889_metawrap_bin.7))+xlab("SRR14722130")+ylab("T1889")

g2<-ggplot(genome_coverage2)+geom_point(aes(x=public_SRR14722130_metawrap_bin.2,y=private_T1916_metawrap_bin.6))+xlab("SRR14722130")+ylab("T1916")

g3<-ggplot(genome_coverage2)+geom_point(aes(x=private_T1889_metawrap_bin.7,y=private_T1916_metawrap_bin.6))+ylab("T1916")+xlab("T1889")


g1 + g2 + g3


colors<-genome_coverage %>% filter(Genome %in% c("public_SRR14722130_metawrap_bin.2","private_T1916_metawrap_bin.6","private_T1889_metawrap_bin.7"))%>%
  pivot_longer(-Genome,values_to="coverage",names_to="metagenome") %>% 
  mutate(presence=ifelse(coverage>50,1,0)) %>% select(-coverage) %>%
  pivot_wider(metagenome,values_from = presence,names_from=Genome) %>%
  mutate(total=public_SRR14722130_metawrap_bin.2+private_T1916_metawrap_bin.6+private_T1889_metawrap_bin.7) %>%
  mutate(col=ifelse(total==0,"has none", ifelse(total==3,"has all",ifelse(total==2 & public_SRR14722130_metawrap_bin.2==0,"has T1916 and T1889",ifelse(total==2 & private_T1916_metawrap_bin.6==0,"has T1889 and SRR14722130",ifelse(total==2,"has T1916 and SRR14722130",ifelse(public_SRR14722130_metawrap_bin.2==1,"has SRR14722130",ifelse(private_T1916_metawrap_bin.6==1,"has T1916","has T1889")))))))) %>%
  select(metagenome,col)
genome_coverage3<-genome_coverage2 %>% left_join(colors)

fig <- plot_ly(genome_coverage3, x = ~public_SRR14722130_metawrap_bin.2, y = ~private_T1889_metawrap_bin.7, z = ~private_T1916_metawrap_bin.6,color=~col,text=~metagenome)
fig


###what is their coverage ratios when they cooccur?
lichenihab_mtg<-mags_role %>% select(Genome,metagenome,depth_cov) %>% 
   filter(Genome %in% c("public_SRR14722130_metawrap_bin.2","private_T1916_metawrap_bin.6","private_T1889_metawrap_bin.7")) %>%
  pivot_wider(names_from=Genome,values_from = depth_cov) 




