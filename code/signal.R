library(tidyverse)
library(plyr)
library(reshape2)
library(AMR)


fungi_genes <- read.table("analysis/05_MAGs/trees/eukaryotes/Fungi/signal_test/fungi_genewise_lnL.txt", check.names=FALSE)
colnames(fungi_genes) <- fungi_genes[1,]
fungi_genes <-fungi_genes[-1,]
count_fungi_genes <- count(fungi_genes$tree_supported) 
colnames(count_fungi_genes)[1] <- "Method"
colnames(count_fungi_genes)[2] <- "n"
count_fungi_genes$Signal <- "Genes"


## Loads site files and set tables 

fungi_site <- read.table("analysis/05_MAGs/trees/eukaryotes/Fungi/signal_test/fungi_sitewise_statictics.table", check.names=FALSE)
colnames(fungi_site)[2] <- "Signal"
colnames(fungi_site)[5] <- "tr1"
colnames(fungi_site)[6] <- "tr2"
fungi_site_m <- melt(fungi_site, id.vars = "Signal", variable.name = "Method", measure.vars = c("tr1", "tr2"))
fungi_site_signal <- fungi_site_m %>% group_by(Signal, Method) %>% filter(Signal != "all") %>% tally(value, name = "n")
fungi_site_all <- fungi_site_m %>% group_by(Signal, Method) %>% filter(Signal == "all") %>% tally(value, name = "n")

fungi_all <- merge(x=fungi_site_all,y=count_fungi_genes, all=TRUE)

fungi_all$Method <- str_replace(fungi_all$Method, "tr1", "Coalescence")
fungi_all$Method <- str_replace(fungi_all$Method, "tr2", "Concatenated")
fungi_all$Signal <- str_replace(fungi_all$Signal, "all", "Sites")
  
## g-test h0: ratio 1:1 
g_fungi_sites<- g.test(fungi_site_all$n)
g_fungi_genes <- g.test(count_fungi_genes$n)


## Barplots

fungi_colors <- c("#4E79A7","#EDC948")
fungi_colors  <- setNames(fungi_colors, levels(fungi_all$Method))

fungi_all %>% arrange(n) %>%
  mutate(Signal = factor(Signal, levels=c("Genes", "Sites")))  %>%
  ggplot(aes(y=n, x=Signal, fill= Method)) + 
  geom_bar(position = "fill", stat="identity") +
  geom_text(aes(label = n),
            stat = "identity",
            colour = "black",
            angle = 90,
            position = position_fill(vjust = 0.5))+
  ylab("Signal Ratio") +
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border =  element_rect(colour = "black", fill=NA),
        panel.background = element_blank()) +
  scale_fill_manual(values = fungi_colors) 

ggsave("analysis/05_MAGs/trees/eukaryotes/Fungi/signal_test/fungi_signal.pdf", width=3, height=4, limitsize = FALSE)


