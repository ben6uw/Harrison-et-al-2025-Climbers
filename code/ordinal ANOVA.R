install.packages('leaflet')

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggpubr)
library(smatr)
library(pheatmap)
library(RColorBrewer)
library(Hmisc)
library(ordinal)
library(igraph)
library(leaflet)


rm(list=ls())
setwd('~/Google Drive/My Drive/Documents/Claudia Sun')

# skip to line ~90+ to read the results of ordinal analysis

load('LCMS_and_biological_data')
# merge LCMS data with biological variables
d[1:7,1:9]


## make a matrix of residuals (mz ~ climbing x age) to call from:
residmat <- matrix(nc=length(mzs), nr=nrow(d))

for(i in 1:length(mzs)){
  residmat[ ,i] <- residuals(lm(d[ ,mzs[i]]~d$Climbing * d$SepDate))} # residuals of main effects of climbing fraction and age (sep date).
colnames(residmat) <- mzs

dat <- data.frame('Climbing'=d$Climbing, residmat)

# an example, not necessarily a 'hit' pair
summary(lm(as.numeric(Climbing) ~ Methionine * AMP, dat)) 

ggplot(dat, aes(AMP, Methionine, color=Climbing))+
  geom_point()+
  theme_bw(base_size = 12)+
  geom_smooth(method='lm', alpha=0.1)+
  facet_wrap(~Climbing, nrow=1, scales='free')


clm(Climbing ~ AMP, data=dat) # cumuative link model, by default y is ordinal and fits the proportional log odds effect of predictors (x) on y
summary(clm(Climbing ~ AMP * Methionine, data=dat))


mzs <- colnames(dat)[-c(1)] # mz names reextracted, to add-back the leading X added

pmat <- matrix(nr=length(mzs), nc=length(mzs))

# this compares a model of non-interaction to interaction on the odds of climbing group: 
# (non-interction) climbing ~  mz1 + mz2
# (interaction) climbing ~  mz1 x mz2, and gives the ANOVA P for the significance of the added interaction term 
for(k in 1:length(mzs)){
  for(i in 1:length(mzs)) {
  a <- anova(clm(Climbing ~ dat[ ,mzs[k]] + dat[ ,mzs[i]], data=dat), clm(Climbing ~ dat[ ,mzs[k]] * dat[ ,mzs[i]], data=dat))
  pmat[k, i] <- a$`Pr(>Chisq)`[2] }}

rownames(pmat) <- mzs
colnames(pmat) <- mzs

save(pmat, dat, d, mzs, residmat, file='mz_mz network/ordinalANOVAresult')



rm(list=ls())
#####################################################################################
# simple summary of ordinal anova result
#####################################################################################
load('mz_mz network/ordinalANOVAresult')

diag(pmat) <- NA
hist(pmat)
table(pmat <=0.05)
pheatmap(pmat, fontsize=5, main='P values')

fdr <- apply(pmat, 2, function(x) p.adjust(x, 'fdr')) # metabolite-wise FDR, correcting for 159 tests for each metabolite

fdr[lower.tri(fdr)] <- NA # this avoids expansion of redundant y~x and x~y by the following line
long <- rstatix::cor_gather(fdr)
head(long)
colnames(long)[3] <- 'FDR'
head(long)
table(is.na(long$FDR))
table(long$FDR <= 0.05) 
n.hits <- sum(long$FDR <= 0.05) # n hits 


rm(list=ls())
#####################################################################################
# permutation testing: are their more edges than expected by chance
#####################################################################################
load('mz_mz network/ordinalANOVAresult')

# test the relevance of climbing fraction on metabolome covariation:
# permutation testing;
# permute samples only, this will:
# preserve within-sample covariance
# shuffle the samples into new climbing fractions

perms <- list()

permat <- matrix(nr=length(mzs), nc=length(mzs))

for(j in 1:100){
  shuf <- sample(nrow(dat), replace=F)
  tmp <- dat[shuf, ]
  for(k in 1:length(mzs)){
  for(i in 1:length(mzs)) {
    a <- anova(clm(dat$Climbing  ~ tmp[ ,mzs[k]] +  tmp[ ,mzs[i]], data=tmp), clm(dat$Climbing  ~ tmp [ ,mzs[k]] *  tmp[ ,mzs[i]], data=tmp))
    permat[k, i] <- a$`Pr(>Chisq)`[2] }}
rownames(permat) <- mzs
colnames(permat) <- mzs
perms[[j]] <- permat }


save(perms, n.hits, d, dat, pmat, residmat, mzs, file='mz_mz network/ordinalANOVA_PermutationTest')
#####################################################################################


rm(list=ls())
#####################################################################################
load('mz_mz network/ordinalANOVA_PermutationTest')

long <- list()

for(j in 1:length(perms)){
fdr <- apply(perms[[j]], 2, function(x) p.adjust(x, 'fdr')) # metabolite-wise FDR, correcting for 159 tests for each metabolite
fdr[lower.tri(fdr)] <- NA # this avoids expansion of redundant y~x and x~y by the following line
long[[j]] <- rstatix::cor_gather(fdr)
colnames(long[[j]])[3] <- paste0('perm_', j, '_FDR')  }

pout <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2), long)

head(pout)
colSums(is.na(pout))
permpairs <- colSums(pout[ ,-c(1:2)] <=0.05)

par(mfrow=c(1,1))

bellringer <- max(c(n.hits, permpairs))
mean(permpairs)
range(permpairs)
hist(permpairs, border = 0, col='grey40', 20, xlim=c(0, (1.15 * bellringer)), las=1, xlab='metabolite pairs (FDR 5%)', main='permutation')
abline(v=n.hits, col=2)
sum(permpairs>=103)/length(permpairs) # empirical P (given 137 pairs in the real data)



rm(list=ls())
############################
# the analysis above only identified metabolites who's covariance depends on ordinal climbing fraction to assess covariance of metabolite pairs by climbing fraction, use major axis regression by climbing fraction:
load('mz_mz network/ordinalANOVAresult')
load('LCMS_and_biological_data') 


# [RECOMMNEND: read this]: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00153.x
pees <- matrix(nr=length(mzs), nc=5) # records P for the mz pairs within each group, could be used later to assess the significance of a given pair in a given climbing group
slopes <- matrix(nr=length(mzs), nc=5) #records the slopes for the mz pairs within each group

# these lists will store slopes and P for each pair within each group, might be useful for a heatmaps perhaps 
Slopelist <- list()
Plist <- list() 

# in the loop below, k = mzY and i = mzX
for(k in 1:length(mzs)) {
  mzY <- mzs[k]
  for (i in 1:length(mzs)) {
    mzX <- mzs[i]
    
    if(mzY == mzX) {pees[i, ]=NA}
    else{m <- sma(residmat[ ,mzY] ~ residmat[ ,mzX] * d$Climbing) # vger
    pees[i, ] <- unlist(m$pval)
    slopes[i, ] <- unlist(m$coef)[grepl(')2', names(unlist(m$coef)))]  }} 
  
  slopes[k, ] =NA # fix for error in loop that propagates stats into situations where mzY = mzX
  
  Plist[[k]] <- pees
  rownames(Plist[[k]]) <- mzs
  rownames(slopes) <- mzs
  Slopelist[[k]] <- slopes
}

names(Plist) <- mzs
names(Slopelist) <- mzs
save(file='sma_results', Plist, Slopelist, d, residmat, mzs)



rm(list=ls())
####################################################

load('mz_mz network/ordinalANOVAresult')
load('mz_mz network/sma_results')
load('LCMS_and_biological_data')

colnames(pmat) <- names(Slopelist)
rownames(pmat) <- names(Slopelist)

diag(pmat) <- NA
hist(pmat)
table(pmat <=0.05)
pheatmap(pmat, fontsize=5, main='P values')

fdr <- apply(pmat, 2, function(x) p.adjust(x, 'fdr')) # metabolite-wise FDR, correcting for 159 tests for each metabolite
fdr[lower.tri(fdr)] <- NA # this avoids expansion of redundant y~x and x~y by the following line
long <- rstatix::cor_gather(fdr)
head(long)
colnames(long)[3] <- 'FDR'
head(long)
table(is.na(long$FDR))
table(long$FDR <= 0.05) # n hits

long[order(long$FDR), ][1:10, ]
summary(sma(X5..Methylthioadenosine ~ Taurine * Climbing, dat))

which(mzs=='Taurine')
Slopelist[[which(mzs=='Taurine')]]["5'-Methylthioadenosine", ]
Plist[[which(mzs=='Taurine')]]["5'-Methylthioadenosine", ]

ggplot(dat, aes(Taurine, X5..Methylthioadenosine, color=Climbing))+
  geom_point()+
  ylab("5'-Methylthioadenosine")+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))+
  geom_smooth(method='lm', alpha=0.1, size=0.7)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  facet_wrap(~Climbing, scales='free', ncol=1)



##################################
### climbing-dependent mz pairs 
##################################
# these will be the edges in the network of each climbign group

mzpairs <- long[long$FDR<=0.05, ]

g <- graph.data.frame(mzpairs[ ,1:2])
plot(g, vertex.size=4, edge.arrow.size=0, vertex.label.cex=0.8, vertex.label.family="Helvetica", vertex.label.dist=1, vertex.label.degree=4, main='climbing group-dependent network')

# to identify edges that are present in top (top and mid-top), and/or bottom (mid-bottom, or bottom), I could use the major axis regression models by climbing group:
load('mz_mz network/sma_results')

maPs <- lapply(Plist, function(x) {rownames(x) <- colnames(pmat); x})
names(maPs) <- mzs

head(mzpairs)
colnames(mzpairs)[3] <- 'ordinalFDR'

# reduce to mzs with sig ordinal Ps
maPs <- maPs[mzpairs$var1] 
maPs

tmp <- maPs
topPs <- numeric()
bottomPs <-  numeric()

for(i in 1:length(maPs)) {
  topPs[i] <- tmp[[i]][mzpairs$var2[i], 1] 
  bottomPs[i] <- tmp[[i]][mzpairs$var2[i], 4]  }

table(topPs<=0.05)
table(bottomPs<=0.05)

ggplot(dat, aes(X6.Methyladenosine, Orotate, color=Climbing))+
  geom_point()+
  theme_bw(base_size = 12)+
  geom_smooth(method='lm', alpha=0.1)+
  facet_wrap(~Climbing, nrow=1, scales='free')

mzpairs$topP <- topPs
mzpairs$bottomP <- bottomPs

top.network  <- graph.data.frame(mzpairs[mzpairs$topP<=0.05 ,1:2])
bottom.network  <- graph.data.frame(mzpairs[mzpairs$bottomP<=0.05 ,1:2])

degree(top.network)[rev(order(degree(top.network)))]
degree(bottom.network)[rev(order(degree(bottom.network)))]

length(V(bottom.network)) # there are 49 pairs in the bottom network
length(V(top.network)) # there are 28 pairs in the bottom network

## use permutation testing to evaluate the significance of the 'bottom' or 'top' metabolite pairs:



#######################################################################################
# compare nodes among top and bottom networks, as well as the main-effect metabolites
#######################################################################################
load('ordinal_main_effect_analysis') # hits from this file are the main-effect mzs from ordinal reg. of climbing on mz
bottom.nodes <- names(V(bottom.network))
top.nodes <- names(V(top.network))

table(hits %in% bottom.nodes)
table(hits %in% top.nodes)
table(bottom.nodes %in% top.nodes)
table(top.nodes %in% bottom.nodes)

# test intersections
fisher.test(table(FDRs[ ,1] <=0.05, mzs %in% bottom.nodes))
fisher.test(table(FDRs[ ,1] <=0.05, mzs %in% top.nodes))

par(mfrow=c(1, 2))
plot(top.network, vertex.size=4, edge.arrow.size=0, vertex.label.cex=0.5, vertex.label.family="Helvetica", vertex.label.dist=1, vertex.label.degree=rep(c(2,4), length(V(top.network))/2), main='top network', vertex.color='blue', layout=layout_on_grid(top.network))

plot(bottom.network, vertex.size=4, edge.arrow.size=0, vertex.label.cex=0.5, vertex.label.family="Helvetica", vertex.label.dist=1, vertex.label.degree=rep(c(2,4), length(V(bottom.network))/2), main='bottom network', vertex.color='red', layout=layout_on_grid(bottom.network))


save(file='mz_mz network/network_data', d, dat, bottom.network, top.network, mzpairs)
###########################################################################################
rm(list=ls())

# make file for cytoscape:
load('mz_mz network/network_data')
# run 'igraphToSif.R' file in /R.functions to get:
source('~My Drive/Documents/R.functions/igraphToSif.R')

igraphToSif(bottom.network, "mz_mz network/bottom.network.sif", "bottom.network")
igraphToSif(top.network, "mz_mz network/top.network.sif", "top.network")



rm(list=ls())
###########################################################################################
# heatmap 
###########################################################################################
load('mz_mz network/network_data')
load('mz_mz network/sma_results')

head(mzpairs)
names(Plist) <- names(Slopelist)

# reduce to mzs with sig ordinal Ps
Slopelist <- Slopelist[mzpairs$var1] 
Slopelist[[1]][1:5, ]
Plist <- Plist[mzpairs$var1] 
Plist[[1]][1:5, ]

tmp <- Slopelist
topSlope <- numeric()
bottomSlope <-  numeric()
allSlopes <- matrix(nr=nrow(mzpairs), nc=4)

for(i in 1:nrow(mzpairs)) {   
  allSlopes[i, ] <- tmp[[i]][mzpairs$var2[i], ]
  topSlope[i] <- tmp[[i]][mzpairs$var2[i], 1] 
  bottomSlope[i] <- tmp[[i]][mzpairs$var2[i], 4]  }

colnames(allSlopes) <- levels(d$Climbing)
rownames(allSlopes) <- paste(mzpairs$var1, '_', mzpairs$var2)


allPs <- matrix(nr=nrow(mzpairs), nc=4)

for(i in 1:nrow(mzpairs)) {   
allPs[i, ] <- Plist[[i]][mzpairs$var2[i], ] }

colnames(allPs) <- levels(d$Climbing)
rownames(allPs) <- paste(mzpairs$var1, '_', mzpairs$var2)


logPs <- -log(allPs, 10)
signPs <- logPs * ifelse(allSlopes>0, 1, -1) # signed -log10P
pheatmap(signPs, fontsize=6, cluster_cols = F, border=0)

my.breaks <- print(seq(-5, 5, by=0.1) )
my.colors <- c(colorRampPalette(colors = c("orange", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "purple"))(length(my.breaks)/2))

pheatmap(signPs, fontsize=6, cluster_cols = F, border=0, color = my.colors, breaks = my.breaks)

rownames(mzpairs) <- paste(mzpairs$var1, '_', mzpairs$var2)
rownames(mzpairs)

top.pairs <- rownames(mzpairs)[mzpairs$topP<=0.05]
bottom.pairs <- rownames(mzpairs)[mzpairs$bottomP<=0.05]

## annotations for heatmap (bottom and top pairs):
top <- as.factor(ifelse(rownames(signPs) %in% top.pairs, 'top', '-'))
bottom <- as.factor(ifelse(rownames(signPs) %in% bottom.pairs, 'bottom', '-'))

annot <- data.frame(top, bottom)
rownames(annot) <- rownames(signPs)
head(annot)

ann_colors =list(top=c(top='red', '-'='grey90'), bottom=c(bottom='#5E3C99', '-'='grey90'))

pheatmap(signPs, fontsize=6, cluster_cols = F, border=0.5, border_color="grey", color = my.colors, breaks = my.breaks, treeheight_row = 0, annotation_row=annot, annotation_colors = ann_colors)
###########

# make split heatmaps, one for top, one for bottom:
topRange <- print(range(signPs[top.pairs, ]))
top.breaks <- seq(topRange[1], topRange[2], by=0.1)
top.colors <- c(colorRampPalette(colors = c("orange", "white"))(length(top.breaks)/2), colorRampPalette(colors = c("white", "purple"))(length(top.breaks)/2))

p1 <- pheatmap(signPs[top.pairs, ], fontsize=6, cluster_cols = F, border=0, color = top.colors, breaks = top.breaks, cellheight = 6)

bottomRange <- range(signPs[bottom.pairs, ])
bottom.breaks <- print(seq(bottomRange[1], bottomRange[2], by=0.1) )
bottom.colors <- c(colorRampPalette(colors = c("orange", "white"))(length(bottom.breaks)/2), colorRampPalette(colors = c("white", "purple"))(length(bottom.breaks)/2))

p2 <- pheatmap(signPs[bottom.pairs, ], fontsize=6, cluster_cols = F, border=0, color = bottom.colors, breaks = bottom.breaks, cellheight = 6)

grid.arrange(grobs = list(p1[[4]], p2[[4]]), ncol=2)
###########



