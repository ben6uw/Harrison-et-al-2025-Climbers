library(igraph)
library(ggraph)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)

rm(list=ls())

dir('~/')
setwd("~/Google Drive/My Drive/Documents/Claudia Sun")

load('mz_mz network/network_data')
mzs <- colnames(d)[7:166]

Climbing.colors <- c("red", "#FDB863", "grey", 'purple', "#5E3C99")
names(Climbing.colors) <- levels(dat$Climbing)
Climbing.colors


bottom.network <- as.undirected(bottom.network)
top.network <- as.undirected(top.network)

table(ifelse(mzpairs$topP <= 0.05, 'top', '-'), ifelse(mzpairs$bottomP <= 0.05, 'bottom', '-')) # could use for venn diagram

length(E(bottom.network))
length(E(top.network))

edge_density(bottom.network)
edge_density(top.network)

degree(bottom.network)[rev(order(degree(bottom.network)))][1:5]
degree(top.network)[rev(order(degree(top.network)))][1:5]

mean_distance(top.network)
mean_distance(bottom.network)

par(mfrow=c(1,1))
plot(bottom.network, vertex.size=4, edge.arrow.size=0, vertex.label.cex=0.5, vertex.label.family="Helvetica", vertex.label.dist=1, vertex.label.degree=rep(c(2,4), length(V(bottom.network))/2), main='bottom network', vertex.color='red', layout=layout_on_grid(bottom.network))

hub_score(top.network, scale=T)$vector # 'hub score' see ?hub.score
hub_score(bottom.network, scale=T)$vector # 'hub score' see ?hub.score

# function for positioning labels outside of the circle:
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (5 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x))) }

top.lab.locs <- radian.rescale(x=1:length(V(top.network)), direction=-1, start=0)
top.hub.scores <- hub_score(top.network, scale=T)$vector

V(top.network)$name <- tolower(V(top.network)$name) # fix up mz names for visualization
V(top.network)$name
V(top.network)$name[V(top.network)$name=='fad'] <- 'FAD'
V(top.network)$name[V(top.network)$name=='ribose-5-p'] <- 'ribose-5-P'
V(top.network)$name[V(top.network)$name=='cysteinyl-glycine (cys-gly)'] <- 'cysteinyl-glycine'
V(top.network)$name[V(top.network)$name=='3-hydroxypropionic acid'] <- '3-hydroxypropioate'
V(top.network)$name[V(top.network)$name=='linoleic acid'] <- 'linolate'
V(top.network)$name[V(top.network)$name=='glutaric acid'] <- 'glutarate'
V(top.network)$name[V(top.network)$name=='citraconic acid'] <- 'citraconate'
V(top.network)$name[V(top.network)$name=='palmitic acid'] <- 'palmitate'
V(top.network)$name[V(top.network)$name=='udp-glcnac'] <- 'UDP-GlcNAc'


plot(top.network, vertex.size=0.005, edge.arrow.size=0, vertex.label.cex=0.6+(top.hub.scores/1.5), vertex.label.family="Helvetica", layout=layout_in_circle (top.network), vertex.label.dist=1.6, vertex.label.degree=top.lab.locs, vertex.color=Climbing.colors[1])


bottom.hub.scores <- hub_score(bottom.network, scale=T)$vector
bottom.lab.locs <- radian.rescale(x=1:length(V(bottom.network)), direction=-1, start=0)

V(bottom.network)$name <- tolower(V(bottom.network)$name) # fix up mz names for visualization
V(bottom.network)$name

V(bottom.network)$name[V(bottom.network)$name=='fad'] <- 'FAD'
V(bottom.network)$name[V(bottom.network)$name=='adp'] <- 'ADP'
V(bottom.network)$name[V(bottom.network)$name=='nad'] <- 'NAD'
V(bottom.network)$name[V(bottom.network)$name=='ribose-5-p'] <- 'ribose-5-P'
V(bottom.network)$name[V(bottom.network)$name=='cysteinyl-glycine (cys-gly)'] <- 'cysteinyl-glycine'
V(bottom.network)$name[V(bottom.network)$name=='3-hydroxypropionic acid'] <- '3-hydroxypropioate'
V(bottom.network)$name[V(bottom.network)$name=='linoleic acid'] <- 'linolate'
V(bottom.network)$name[V(bottom.network)$name=='glutaric acid'] <- 'glutarate'
V(bottom.network)$name[V(bottom.network)$name=='citraconic acid'] <- 'citraconate'
V(bottom.network)$name[V(bottom.network)$name=='palmitic acid'] <- 'palmitate'
V(bottom.network)$name[V(bottom.network)$name=='udp-glcnac'] <- 'UDP-GlcNAc'
V(bottom.network)$name[V(bottom.network)$name=='dimethylarginine (a/sdma)'] <- 'dimethylarginine'
V(bottom.network)$name[V(bottom.network)$name=='2-aminoisobutyric acid'] <- '2-aminoisobutyrate'
V(bottom.network)$name[V(bottom.network)$name=='n-acetylglycine'] <- 'n-Ac-glycine'
V(bottom.network)$name[V(bottom.network)$name=='n-ac-alanine'] <- 'n-Ac-alanine'
V(bottom.network)$name[V(bottom.network)$name=='n-ac-arginine'] <- 'n-Ac-arginine'
V(bottom.network)$name[V(bottom.network)$name=='glycerol-3-p'] <- 'glycerol-3-P'
V(bottom.network)$name[V(bottom.network)$name=='dl dopa'] <- 'DL-dopa'
V(bottom.network)$name[V(bottom.network)$name=='phenylacetic acid'] <- 'phenylacetate'
V(bottom.network)$name[V(bottom.network)$name=='l-kynurenine'] <- 'L-kynurenine'


bottom.lab.locs <- radian.rescale(x=1:length(V(bottom.network)), direction=-1, start=0)

plot(bottom.network, vertex.size=0.005, edge.arrow.size=0, vertex.label.cex=0.6+hub_score(bottom.network, scale=F)$vector, vertex.label.family="Helvetica", layout=layout_in_circle (bottom.network), vertex.label.dist=1.6, vertex.label.degree=bottom.lab.locs, vertex.color=Climbing.colors[5])

plot(top.network, vertex.size=0.005, edge.arrow.size=0, vertex.label.cex=0.6+(top.hub.scores/1.5), vertex.label.family="Helvetica", layout=layout_in_circle (top.network), vertex.label.dist=1.6, vertex.label.degree=top.lab.locs, vertex.color=Climbing.colors[1])


#########################################
# follow these steps to export to Cytoscape
BiocManager::install('RCy3')

require(RCy3)
# 1. Download the latest Cytoscape from http://www.cytoscape.org/download.php
# 2. Complete installation wizard
# 3. Launch Cytoscape
cytoscapePing() # confirm that this line gives: "You are connected to Cytoscape!"

createNetworkFromIgraph(bottom.network, title = "bottom.network", collection = "My Igraph Network Collection") # this should open in Cytoscape automatically

createNetworkFromIgraph(top.network, title = "top.network", collection = "My Igraph Network Collection") # this should open in Cytoscape automatically

#[NOTE]: in Cytoscape: File>Import>Table from File  to import the table from generateResultsTable(). This will allow you to import the node types (i.e. pathways, reactions, enzymes, etc.) into the Node Table in Cytopscape.
#########################################


# plot node degree by metabolite:
degb <- as.data.frame(degree(bottom.network, v = V(bottom.network), mode = c("all"),
       loops = TRUE, normalized = FALSE))

colnames(degb) <- 'degree'
degb$mz <- rownames(degb)
degb
degb$mz <- factor(degb$mz, levels=degb$mz[rev(order(degb$degree))])

bdd <- ggplot(degb, aes(y=degree, x=mz))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('node degree')+
  xlab('')


degt <- as.data.frame(degree(top.network, v = V(top.network), mode = c("all"),
                             loops = TRUE, normalized = FALSE))

colnames(degt) <- 'degree'
degt$mz <- rownames(degt)
degt
degt$mz <- factor(degt$mz, levels=degt$mz[rev(order(degt$degree))])

tdd <- ggplot(degt, aes(y=degree, x=mz))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('node degree')+
  xlab('')

ggarrange(tdd, bdd, ncol=1)



degt$grp <- 'top'
degb$grp <- 'bottom'
x <- rbind(degb, degt)

head(x)

x <- pivot_wider(x, values_from = degree, names_from = grp)

plot(x$bottom, x$top, pch=19, col=rgb(0, 0, 0, 0.5)) # no evidence of a relationship btw degrees of a particular mz in bottom and top network

             