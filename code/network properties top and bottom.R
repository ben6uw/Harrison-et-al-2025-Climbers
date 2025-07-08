library(igraph)
library(ggraph)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)

rm(list=ls())

dir('~/')
setwd("~/Google Drive/My Drive/Documents/Claudia Sun")

load('mz_mz network/network_data') # file generated in 'ordinal ANOVA.R'
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
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x))) }

top.lab.locs <- radian.rescale(x=1:length(V(top.network)), direction=-1, start=0)
top.hub.scores <- hub_score(top.network, scale=T)$vector

plot(top.network, vertex.size=5*top.hub.scores, edge.arrow.size=0, vertex.label.cex=0.4+top.hub.scores, vertex.label.family="Helvetica", layout=layout_in_circle (top.network), vertex.label.dist=1.6, vertex.label.degree=top.lab.locs, vertex.color=Climbing.colors[1])

sort(top.hub.scores, decreasing =T) # scores are scaled


bottom.hub.scores <- hub_score(bottom.network, scale=T)$vector
bottom.lab.locs <- radian.rescale(x=1:length(V(bottom.network)), direction=-1, start=0)

plot(bottom.network, vertex.size=10*(hub_score(bottom.network, scale=T)$vector), edge.arrow.size=0, vertex.label.cex=0.4+hub_score(bottom.network, scale=F)$vector, vertex.label.family="Helvetica", layout=layout_in_circle (bottom.network), vertex.label.dist=1.6, vertex.label.degree=bottom.lab.locs, vertex.color=Climbing.colors[5])

sort(bottom.hub.scores, decreasing =T) # scores are scaled

  