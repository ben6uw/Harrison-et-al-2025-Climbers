library(igraph)
library(ggplot2)
library(tidyverse)
library(data.table)
library(limma)

# mediator analysis for 'bottom' and 'top' networks:
rm(list=ls())

setwd('~/Google Drive/My Drive/Documents/Claudia Sun')

load('LCMS_and_biological_data') # save this version of 'd'
table(d$Climbing)
d2 <- d
load('mz_mz network/network_data') # file generated in 'ordinal ANOVA.R'
# reloads a smaller version of 'd' without the middle climbing group.
rm(d)
d <- d2
load('path.analysis.2') # loads allpair, other pairs

save(allpair, bottom.network, d, d2, graph, meds, mz.info, mzpairs, path.nodes, top.network, kegg, keggmzs, mzs, sepDays, file='dataForMediationAnalysis')



head(allpair)
mean(allpair$total.nodes)
mean(allpair$unique.nodes)

p1 <- ggplot(allpair, aes(x=total.nodes)) + geom_histogram() + scale_x_log10(breaks=c(1,2,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('total KEGG nodes among all shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')

p2 <- ggplot(allpair, aes(x=unique.nodes)) + geom_histogram() + scale_x_log10(breaks=c(1,2,5,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('unique KEGG nodes among all shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')

p3 <- ggplot(allpair, aes(x=n.paths)) + geom_histogram() + scale_x_log10(breaks=c(1,2,5,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('number of shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')


ggarrange(p3, p1, p2, ncol=1)


table(allpair$mz_mz_fdr<=0.05)
allpair$sig <- as.factor(ifelse(allpair$mz_mz_fdr<=0.01, 'FDR<0.01', 'ns'))

ggplot(allpair, aes(x=total.nodes, fill=sig)) + geom_histogram() + scale_x_log10(breaks=c(1,2,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('total KEGG nodes among all shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')+
  scale_fill_manual(values=c('purple', 'grey'))+
  facet_wrap(~sig, scales='free', ncol=1)

ggplot(allpair, aes(x=unique.nodes, fill=sig)) + geom_histogram() + scale_x_log10(breaks=c(1,2,5,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('unique KEGG nodes among all shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')+
  scale_fill_manual(values=c('purple', 'grey'))+
  facet_wrap(~sig, scales='free', ncol=1)


ggplot(allpair, aes(x=n.paths, fill=sig)) + geom_histogram() + scale_x_log10(breaks=c(1,2,5,10,100, 1000, 10000, 50000))+
  theme_bw(base_size = 12)+
  xlab('number of shortest paths between pair')+
  ylab('metabolite pairs (n=8,778)')+
  scale_fill_manual(values=c('purple', 'grey'))+
  facet_wrap(~sig, scales='free', ncol=1)







rm(list=ls())
#######################################################################
load('dataForMediationAnalysis')



## summary of networks:
#######################################################################
# bottom.netwrok = all pairs of targeted metabolites correlated at P<0.05 in the bottom fraction flies
# top.network = all pairs of targeted metabolites correlated at P<0.05 in the top fraction flies
# all.pair = all pairs of targeted metabolites that were present in the KEGG graph (some targeted metaboites have no KEGG ID, others are not associated with any nodes on the KEGG graph)

sqr <- length(mzs)^2
sqr/2 - 0.5*length(mzs) # number of pairs of mzs that are possible in the 160-metabolites measured on the targeted LCMS = 12,720    NOTE: the KEGG network only has 133 of the mzs, which can then only make 8,778 possible pairs

# convert networks to undirected
bottom.network <- as.undirected(bottom.network)
bottompairs <- as.data.frame(as_edgelist(bottom.network))
nrow(bottompairs) # 55 pairs in the bottom network

# join the bottom network with the pairs of metabolites in KEGG:
colnames(bottompairs) <- c('row', 'col')
bottom1 <- plyr::match_df(allpair, bottompairs)
colnames(bottompairs) <- c('col', 'row')
bottom2 <- plyr::match_df(allpair, bottompairs)
allpair$bottom <- ifelse(rownames(allpair) %in% c(rownames(bottom1), rownames(bottom2)), 'bottom', 'bottombkd')


## summary from the bottom network:
table(mzs %in% c(bottompairs$col, bottompairs$row)) # among 49 targeted metabolites among the pairs in the bottom network, of these:
length(unique(c(allpair$kegg1[allpair$bottom=='bottom'], allpair$kegg2[allpair$bottom=='bottom']))) # 31 metabolites in the bottom network are on the KEGG graph, and they form:
nrow(bottompairs) # 55 pairs on the bottom covariance network, and of these:
table(allpair$bottom) # 29 pairs are on the KEGG graph
head(allpair)

# 31 bottom metabolites are on the KEGG network, and they make 29 pairs
allpair[allpair$bottom == 'bottom', c('kegg1',  'kegg2')]


#######################################################################

top.network <- as.undirected(top.network)
toppairs <- as.data.frame(as_edgelist(top.network))
nrow(toppairs) 

colnames(toppairs) <- c('row', 'col')
top1 <- plyr::match_df(allpair, toppairs)
colnames(toppairs) <- c('col', 'row')
top2 <- plyr::match_df(allpair, toppairs)
allpair$top <- ifelse(rownames(allpair) %in% c(rownames(top1), rownames(top2)), 'top', 'topbkd')

## summary from the top network:
table(mzs %in% c(toppairs$col, toppairs$row)) # among 28 targeted metabolites among the pairs in the top network
length(unique(c(allpair$kegg1[allpair$top=='top'], allpair$kegg2[allpair$top=='top']))) # 22 are on the KEGG network, and they form:
table(allpair$top) # 21 pairs are on the KEGG graph

#######################################################################
#######################################################################





rm(list=ls())

################################################################################################
# 'mediator analysis', looking for enrichment of kegg entities (ie. pahtways or enzymes, etc.) among nodes in shortest paths connecting metabolites of interest:
################################################################################################
load('data.for.mediation.analysis')

##############################################################
# test for mediating PATHWAYS
kegg.pathways <- names(V(graph))[grepl('dme', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.pathways)) {
  if (length(unique(path.nodes %like% kegg.pathways[i]))>1) {
    tab <- table(ifelse(allpair$bottom =='bottom', 'pair', 'not'), ifelse(path.nodes %like% kegg.pathways[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}

table(is.na(p)) # 77 of 133 kegg pathways are part of at least one shortest path
# 56 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.pathways
names(sign) <- kegg.pathways

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.pathways
table(fdr <=0.05)

P.meds <- data.frame('pathway'=kegg.pathways, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- P.meds$pathway[P.meds$moreThanExpected >0 & P.meds$chiFDR <=0.05]
mediators <- print(mediators[!is.na(mediators)])

save(P.meds, file='mediationResults/mediation_bottom_PATHWAYS')
##############################################################


##############################################################
# test for mediating MODULES
kegg.modules <- names(V(graph))[grepl('M0', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.modules)) {
  if (length(unique(path.nodes %like% kegg.modules[i]))>1) {
    tab <- table(ifelse(allpair$bottom =='bottom', 'pair', 'not'), ifelse(path.nodes %like% kegg.modules[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}


table(is.na(p)) # 14 of 174 kegg modules are part of at least one shortest path
# 160 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.modules
names(sign) <- kegg.modules

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.modules
table(fdr <=0.05)

M.meds <- data.frame('pathway'=kegg.modules, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- P.meds$pathway[P.meds$moreThanExpected >0 & P.meds$chiFDR <=0.05]
mediators <- print(mediators[!is.na(mediators)])

save(M.meds, file='mediationResults/mediation_bottom_MODULES')
##############################################################


##############################################################
# test for mediating ENZYMES
kegg.enzymes <- names(V(graph))[grepl('\\.', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.enzymes)) {
  if (length(unique(path.nodes %like% kegg.enzymes[i]))>1) {
    tab <- table(ifelse(allpair$bottom =='bottom', 'pair', 'not'), ifelse(path.nodes %like% kegg.enzymes[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}


table(is.na(p)) # 202 of 768 kegg enzymes are part of at least one shortest path
# 566 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.enzymes
names(sign) <- kegg.enzymes

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.enzymes
table(fdr <=0.05)

E.meds <- data.frame('pathway'=kegg.enzymes, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- P.meds$pathway[P.meds$moreThanExpected >0 & P.meds$chiFDR <=0.05]
mediators <- print(mediators[!is.na(mediators)])

save(E.meds, file='mediationResults/mediation_bottom_ENZYMES')
##############################################################

## compile results
load('mediationResults/mediation_bottom_PATHWAYS')
load('mediationResults/mediation_bottom_MODULES')
load('mediationResults/mediation_bottom_ENZYMES')

mediationResultsList <- list(P.meds, M.meds, E.meds)
names(mediationResultsList) <- c('PATHWAYS', 'MODULES', 'ENZYMES')

mediationResults <- as.data.frame(do.call(rbind, mediationResultsList))
mediationResults$mediatorType <- rownames(mediationResults)
mediationResults$network <- 'bottom'

# extract mediators
meds <- P.meds[P.meds$moreThanExpected>0 & P.meds$chiFDR <=0.05, ]
meds <- meds[rowSums(is.na(meds))==0, ]


a <- limma::getKEGGPathwayNames(species="dme")
head(a)
colnames(a)[1] <- 'pathway'

meds <- merge(meds, a, all.x=T)

meds
write.csv(meds, quote=F, row.names = F, file='bottom.hits.csv')



rm(list=ls()) # now for enrichement among top network: 
################################################################################################
load('data.for.mediation.analysis') 

# mediating pathways?
kegg.pathways <- names(V(graph))[grepl('dme', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.pathways)) {
  if (length(unique(path.nodes %like% kegg.pathways[i]))>1) {
    tab <- table(ifelse(allpair$top =='top', 'pair', 'not'), ifelse(path.nodes %like% kegg.pathways[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}


table(is.na(p)) # 77 of 133 kegg pathways are part of at least one shortest path
# 56 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.pathways
names(sign) <- kegg.pathways

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.pathways
table(fdr <=0.05)

min(p[!is.na(p)])
min(fdr[!is.na(fdr)])

P.meds <- data.frame('pathway'=kegg.pathways, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- P.meds$pathway[P.meds$moreThanExpected >0 & P.meds$chiFDR <=0.05]
mediators <- mediators[!is.na(mediators)]



# mediating modules?
kegg.modules <- names(V(graph))[grepl('M0', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.modules)) {
  if (length(unique(path.nodes %like% kegg.modules[i]))>1) {
    tab <- table(ifelse(allpair$top =='top', 'pair', 'not'), ifelse(path.nodes %like% kegg.modules[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}


table(is.na(p)) # 14 of 174 kegg modules are part of at least one shortest path
# 160 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.modules
names(sign) <- kegg.modules

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.modules
table(fdr <=0.05)

min(p[!is.na(p)])

M.meds <- data.frame('pathway'=kegg.modules, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- M.meds$pathway[M.meds$moreThanExpected >0 & M.meds$chiFDR <=0.05]
mediators <- mediators[!is.na(mediators)]



# mediating enzymes?
kegg.enzymes <- names(V(graph))[grepl('\\.', names(V(graph)))]
p <- numeric()
sign <- numeric()
mediator.tabs <- list()

for(i in 1:length(kegg.enzymes)) {
  if (length(unique(path.nodes %like% kegg.enzymes[i]))>1) {
    tab <- table(ifelse(allpair$top =='top', 'pair', 'not'), ifelse(path.nodes %like% kegg.enzymes[i], 'mediator', 'not'))
    p[i] <- fisher.test(tab)$p.value 
    mediator.tabs[[i]] <- tab 
    cst <- chisq.test(tab)
    sign[i] <- (cst$observed - cst$expected)[2,1]} 
  else { 
    p[i] <- NA 
    sign[i] <- NA
    mediator.tabs[[i]] <- NA
  }}


table(is.na(p)) # 202 of 768 kegg enzymes are part of at least one shortest path
# 566 kegg pathways are not on a shortest path

names(mediator.tabs) <- kegg.enzymes
names(sign) <- kegg.enzymes

fdr <- p.adjust(p, 'fdr')
names(fdr) <- kegg.enzymes
table(fdr <=0.05)

min(p[!is.na(p)])

E.meds <- data.frame('pathway'=kegg.enzymes, 'moreThanExpected'=sign, 'chiP'=p, 'chiFDR'=fdr)
mediators <- E.meds$pathway[E.meds$moreThanExpected >0 & E.meds$chiFDR <=0.05]
mediators <- mediators[!is.na(mediators)]

meds <- P.meds[P.meds$moreThanExpected>0 & P.meds$chiFDR <=0.05, ]
meds <- meds[rowSums(is.na(meds))==0, ]

a <- getKEGGPathwayNames(species="dme")
head(a)
colnames(a)[1] <- 'pathway'

meds <- merge(meds, a, all.x=T)

meds
write.csv(meds, quote=F, row.names = F, file='top.hits.csv')




rm(list=ls()) 
##################################################################
# plotting a sub network to depict the mediation enrichemnent:
##################################################################
load('data.for.mediation.analysis') 
hits <- read.csv('bottom.hits.csv')

hits
hits$pathway

# "dme00430" of the 29 pairs in the bottom network, "dme00430" is in 5 of those paths. The 5 paths occupied by "dme00430" in the bottom network are 5/46, (a little over 10%) of all pairs with short paths containing "dme00430"

# since a hit connects all pairs-of-interest by definition, could make a network to show what it explains:

# load paths
# chop down to only those among bottom pairs
# then chop to those that contain the hit:

head(allpair)

sk1 <- allpair$kegg1[allpair$bottom =='bottom']
sk2 <- allpair$kegg2[allpair$bottom =='bottom']

pathList <- list()
for(i in 1:length(sk1)) {
pathList[[i]] <- all_shortest_paths(graph, sk1[i], sk2[i], weights=NA)  }

mediated.pairs <- unlist(lapply(lapply(pathList, function(x) unique(names(unlist(x$res)))), function(x) "dme00430"%in% x))


subList <- (pathList)[mediated.pairs] # reduce to only pairs with 'hit' among their paths
members <- lapply(subList, function(x) unlist(lapply(x$res, function(x) 'dme00430' %in% names(x)))) # paths mediated in each sublist element
members

subplots <- list()
length(subList) # number of mediated pairs
dfs <- list()

for(k in 1:length(subList) ) { # k will be an mz pair
vertList <- print((subList[[k]]$res)[unlist(members[k])]) # gets all paths that contain the 'hit' from mzpair k

dfList <- list()
for(j in 1:length(vertList)) {
verts <- names(unlist(vertList[j]))
len <- length(verts)
index <- 1:len

tmp <- matrix(nr=(len-1), nc=2)
for(i in 1:(len-1)){
tmp[i, ] <- index[c(i, i+1)] }
dfList[[j]] <- data.frame(verts[tmp[ ,1]], verts[tmp[ ,2]]) }

df <- do.call(rbind, dfList)
df <- df[!duplicated(df), ]

dfs[[k]] <- df
subplots[[k]] <- graph_from_data_frame(df, directed = F) }


lapply(subplots, length) # the size of each igraph subnetwork
plot(subplots[[1]]) # an example
plot(subplots[[5]]) # an example

# plot each subnetwork:
tmpVertNames <- c(hits$pathway, keggmzs)
names(tmpVertNames)[1] <- hits$pathway

par(mfrow=c(3,2))
for(k in 1:length(subplots)){
tmp <- subplots[[k]]
V(tmp)$label [names(V(tmp)) %in% names(tmpVertNames)] <- tmpVertNames[names(V(tmp))[names(V(tmp)) %in% names(tmpVertNames)]]

plot(tmp, vertex.size=hub_score(tmp, scale=T)$vector, edge.arrow.size=0, vertex.label.cex=1, vertex.label.family="Helvetica", layout=layout_with_kk(tmp), vertex.label.dist=1) }
#######################

## for visualization:
# build ONE network that includes all pairs mediated by dme00430:
cbind(keggmzs[sk1][mediated.pairs], keggmzs[sk2][mediated.pairs]) # the 5 pairs of mzs that are mediated by dme00430

df2 <- do.call(rbind, dfs)
df2
df2 <- df2[!duplicated(df2), ] # remove duplicate pairs
sub <- graph_from_data_frame(df2, directed = F) # build a common network - all roads go through dme00430
plot(sub) 

tmpVertNames <- c(hits$pathway, tolower(keggmzs)) # fix up mz names for visualization
tmpVertNames[tmpVertNames =='nad'] <- 'NAD'
names(tmpVertNames)[1] <- hits$pathway

V(sub)$label [names(V(sub)) %in% names(tmpVertNames)] <- tmpVertNames[names(V(sub))[names(V(sub)) %in% names(tmpVertNames)]]

par(mfrow=c(1,1))  

plot(sub, vertex.size=3, edge.arrow.size=0, vertex.label.cex=1, vertex.label.family="Helvetica", layout=layout_with_kk(sub), vertex.label.dist=1, vertex.color = c("white", "tomato", "purple")[1+ names(V(sub)) %in% names(keggmzs) + 2*names(V(sub)) %in% hits$pathway]) 


V(sub)$name
V(sub)$label[!V(sub)$name %in% names(tmpVertNames)] <- V(sub)$name[!V(sub)$name %in% names(tmpVertNames)]
   

plot(sub, vertex.size=4, edge.arrow.size=0, vertex.label.cex=1, vertex.label.family="Helvetica", layout=layout_with_kk(sub), vertex.label.dist=1, vertex.color = c("grey", "tomato", "purple", 'lightblue', 'blue')[1+ names(V(sub)) %in% names(keggmzs) + 2*grepl('dme', names(V(sub))) + 3*grepl('C', names(V(sub))) + 4*grepl('\\.', names(V(sub)))])


otherCompounds[V(sub)$label]

neems <- c("C00080", "C00005", "C00006", "C00007", "C00004", "C00009", "C00001", "C00014", "C00013")
otherCompounds <- c('H+', 'NADPH', 'NADP+', 'oxygen', 'NADH', 'Phosphate', 'water', 'ammonia', 'diphosphate')
names(otherCompounds) <- neems
otherCompounds



rm(list=ls())
##################################################################
# testing whether a metabolite set (ie, all mzs in bottom network), and not necessarily the covariance between those mzs, cold explain node enrichment
##################################################################
load('data.for.mediation.analysis') 

plot(bottom.network, vertex.size=1,  vertex.label.cex=0.2)

r <- bottom.network # for ease of coding, 'r' is the real network
w <- r
V(w)$name=sample(V(r)$name) # simply shuffle node names
par(mfrow=c(1,2))
plot(r, vertex.size=1,  vertex.label.cex=1)
plot(w, vertex.size=1,  vertex.label.cex=1)


l <- as_long_data_frame(r)
head(l)

# join the bottom network with the pairs of metabolites in KEGG:
l <- as_long_data_frame(r)[ ,3:4]
colnames(l) <- c('row', 'col')
l1 <- plyr::match_df(allpair, l)
colnames(l) <- c('col', 'row')
l2 <- plyr::match_df(allpair, l)
allpair$bottom <- ifelse(rownames(allpair) %in% c(rownames(l1), rownames(l2)), 'bottom', 'bottombkd')

table(allpair$bottom)

head(allpair)
x <- allpair
x$dme00430inPath <- path.nodes %like% 'dme00430'

x <- x[x$bottom=='bottom', ] 
table(x$path, x$dme00430inPath)


##############################################################
table(path.nodes %like% 'dme00430') # dme00430 is present among the paths that connect 46 pairs of mzs on KEGG.

# review the result in the real network:
tab <- table(ifelse(allpair$bottom =='bottom', 'pair', 'not'), ifelse(path.nodes %like% 'dme00430', 'mediator', 'not')) # make contingency table, if a pair covaries in the bottom network, how often does it also contain dme00430?

tab
fisher.test(tab) # use Fisher test, rather than Chi Sq, as the expectations can be low-freq

# permutation test
tabList <- list()
mediate <- ifelse(path.nodes %like% 'dme00430', 'mediator', 'not')

nperms <- 10000

for(i in 1:nperms) {
V(w)$name=sample(V(r)$name) # simply shuffle node names
l <- as_long_data_frame(w)[ ,3:4]
colnames(l) <- c('row', 'col')
l1 <- plyr::match_df(allpair, l)
colnames(l) <- c('col', 'row')
l2 <- plyr::match_df(allpair, l)
allpair$perm.bottom <- ifelse(rownames(allpair) %in% c(rownames(l1), rownames(l2)), 'perm.bottom', 'perm.bottombkd')
pair <- ifelse(allpair$perm.bottom =='perm.bottom', 'pair', 'not')
tabList[[i]] <- table(pair, mediate) }

present <- sapply(tabList, function(x) return(x[2,1]))

par(mfrow=c(1,1))
hist(present, 12, col='grey40', border=0, xlab='N paths containing dme00430', main='')
abline(v=5, col='red')
mean(present)
median(present)
range(present)

tab # the real result
s <- sum(present > tab[2, 1]) # number of times that the permuted network had a greater number of mediated pairs than the real netwrok
(s+1)/(nperms+1) # empirical P value

