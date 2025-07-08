
# Analyze metabolite dat from Claudia's experiment

# BiocManager::install('sva')

suppressPackageStartupMessages({
  library(data.table)  
  library(tidyverse)  
  library(stringr)  
  library(RColorBrewer)  
  library(lme4)  
  library(lmerTest)  
  library(sva)
  library(caret)  
  library(ggpubr)  
  library(RColorBrewer)  
  library(pheatmap)  
  library(sjPlot)  
  library(ggforce)
  library(sjmisc)  
  library(pls)  
  library(igraph)  
  library(dendextend)  
  library(ggrepel)  
  library(effects)  
  library(emmeans)  
  library(effectsize)  
  library(sjstats)  
  library(gridExtra)  
  library(knitr)  
  library(rmarkdown)  
  library(MASS)
  library(ordinal)
  library(brant)
})

rm(list=ls())

getwd()
setwd("/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/Claudia Sun")

## proceed to line ~200 if you have already normalized the data

dat <- read.table("Claudias data and code/compound.csv", sep= ",", quote = "\"", header = T, check.names = F)
vials <- read.table("Claudias data and code/Vial-Specific info.csv", sep= ",", header = T)
colnames(vials)[1] <- 'Vial'
samps<- read.csv("Claudias data and code/sample order.csv", header = T)
QC <- read.csv("Claudias data and code/QC.csv", quote = "\"", header = T, check.names = F)
BI <- read.csv("Claudias data and code/Batch info.csv", quote = "\"", header = T, check.names = F)
dat <- merge(dat, QC, by.x= "COMPOUND", by.y= "Current MS Compounds")

# missingness, NAs in the raw dat
dat[1:4,1:4]
dat[dat=="N/A"]=NA
mzdat <- as.matrix(dat[ ,4:ncol(dat)])

colnames(mzdat) # samples 1:72 are biosamples, the remainder are QCs
hist(rowSums(is.na(mzdat[ ,1:75])))
table(rowSums(is.na(mzdat[ ,1:72])) ==0) # 162 complete mzs
table(rowSums(is.na(mzdat[ ,1:72]))) # the next most-complete mz is missign in 6 of 72 samples

# There are two control metabolites are included in those 162 metabolites
dat <- dat[rowSums(is.na(dat[ ,4:75])) == 0, ]
table(grepl('13', dat$COMPOUND)) # 2 spiked isotopes among the metabolites
mzs <- dat$COMPOUND[!grepl('13', dat$COMPOUND)] #160 metabolites
mz.info <- dat[dat$COMPOUND %in% mzs, 1:3] # pull KEGG and HMDB identities

samples <- colnames(dat)[-c(1:3)]
rownames(dat) <- dat$COMPOUND

mzdat <- t(dat[ ,4:ncol(dat)])
mzdat <- mzdat[ ,mzs]
mzdat <- apply(mzdat, 2, as.numeric)
rownames(mzdat) <- samples
mzdat[1:4,1:4]
rowSums(is.na(mzdat[ ]))

#merge with sample order
head(samps)
rownames(samps) <- samps$Sample.ID
d <- cbind(samps, mzdat[samps$Sample.ID, ])
d[1:6,1:10]

#merge with Batch Info
head(BI)
d$Sample.ID[d$Sample.ID %in% BI$Vial]
d$Sample.ID[!d$Sample.ID %in% BI$Vial]
d <- merge(BI, d, by.x= "Vial", by.y="Sample.ID", all.y=T)
d$Batch <- as.factor(d$Batch)
d[1:7,1:5]

# transform 
par(mar=c(5,4,4,2))
par(mfrow=c(2,2))

boxplot(d[ ,mzs], las=1, ylab='raw (m/z)')
boxplot(log(d[ ,mzs]), las=1, ylab='log(m/z)')
d[ ,mzs] <- log(d[ ,mzs])

boxplot(t(d[ ,mzs]), col=ifelse(grepl('QC', d$Vial), 2, 'grey'), las=1, xlab='samples', ylab='log(m/z) data')
sc <- apply(d[ ,mzs], 1, function(x) scale(x, scale=F)) # mean center the data, but do not scale
boxplot(sc, col=ifelse(grepl('QC', d$Vial), 2, 'grey'), las=1, xlab='samples', ylab='mean-centered log(m/z) data')
d[ ,mzs] <- t(sc) # replace data with scaled logmzs

d$sampleType <- d$Vial
d$sampleType[!grepl('QC', d$sampleType )] <- 'fly.sample'
d$sampleType[grepl('S', d$sampleType )] <- 'pooled.sample'
d$sampleType[grepl('I', d$sampleType )] <- 'instrumentQC'
d$sampleType <- as.factor(d$sampleType)
table(d$sampleType)

# plot metabolites by run order
colnames(d)[1:10]
long <- d %>% gather(metabolite, log_mz, -c(Vial, Batch, Sample.order, sampleType))
head(long)

ggplot(subset(long, metabolite %in% mzs[1:20]), aes(y=log_mz, x=Sample.order, color=sampleType))+
  geom_point()+
  theme_bw()+
  facet_wrap(~metabolite)


tmp <- d
d <- d[!grepl('QC', d$Vial ), ] # remove QC samples prior to removing batch effects, and run order effects

par(mar=c(5,4,4,2))
par(mfrow=c(2,2))

plot(prcomp(d[ ,mzs], scale=T)$x, col=d$Batch, pch=16, main='metabolite extraction batch', las=1) # big ol batch effect

comb <- sva::ComBat(t(d[ ,mzs]), batch=d$Batch)
plot(prcomp(t(comb), scale=T)$x, col=d$Batch, pch=16, main='combat vs. batch', las=1)
d[ ,mzs] <- t(comb) # accept batch correction 


# run order effects?
stripchart(Sample.order ~ Batch, d, col=1, pch=1, las=1, ylab='prep batch', main='prep batches across LCMS run')

# add QC data back into data and plot all metabolites by run order
tmp <- rbind(d, tmp[grepl('QC', tmp$Vial ), ])
long <- tmp %>% gather(metabolite, log_mz, -c(Vial, Batch, Sample.order, sampleType))
head(long)

ggplot(subset(long, metabolite %in% mzs[1:20]), aes(y=log_mz, x=Sample.order, color=sampleType))+
  geom_point()+
  theme_bw()+
  facet_wrap(~metabolite)


# analyze run order effects per metabolite (do this on the data WITHOUT the QCs):
p <- numeric()
for (i in 1:length(mzs)) {
  s<- summary(lm(d[ ,mzs[i]] ~ d$Sample.order))
  p[i] <- s$coefficients[2,4] }

FDR <- p.adjust(p, 'BH')  
runO <- as.data.frame(cbind(p, FDR))
rownames(runO) <- mzs
head(runO[order(runO$p), ])
plot(Creatine ~ Sample.order, d, pch=16, main='run order effect')

mois <- rownames(runO)[runO$FDR < 0.05] # some of the most dramatic offenders

par(mfrow=c(5,4))
par(mar=c(2,2,2,2))
for (i in 1:length(mois)) {
  plot(d[ ,mois[i]] ~ d$Sample.order, pch=16, cex=0.5, xlab='run order', ylab='', las=1, main=mois[i])
  abline(lm(d[ ,mois[i]] ~ d$Sample.order)) }

for(i in 1:7) {
plot(prcomp(d[ ,mzs], scale=T)$x[ ,i] ~ d$Sample.order, pch=16, main=paste('PC', i), cex=0.5, col=2)} # run order effects on PC2, PC4

# how to handle run order effects?

# regress-out mz ~ run order from all metabolites
tmp <- d

par(mfrow=c(2,3))
par(mar=c(5,5,4,2))
boxplot(t(d[ ,mzs]), las=1, ylab='batch-corrected log(mz)', xlab='sample', main='pre run order correction', cex.lab=2)
for (i in 1:length(mzs)) {
  tmp[ ,mzs[i]] <- lm(d[ ,mzs[i]] ~ d$Sample.order)$residuals }
boxplot(t(tmp[ ,mzs]), main='mz ~ run order residuals', ylab='~ run order residuals', xlab='sample', las=1, cex.lab=2) 

plot(apply(d[ ,mzs], 1, var), apply(tmp[ ,mzs], 1, var), xlab='sample variance (pre correction)',  ylab='sample variance (post correction)', las=1, pch=16, cex.lab=2)

boxplot(d[ ,mzs], las=1, ylab='batch-corrected log(mz)', xlab='metabolite', main='pre run order correction', cex.lab=2)
boxplot(tmp[ ,mzs], main='mz ~ run order residuals', ylab='~ run order residuals', xlab='metabolite', las=1, cex.lab=2) 
plot(apply(d[ ,mzs], 2, var), apply(tmp[ ,mzs], 2, var), xlab='metabolite variance (pre correction)',  ylab='metabolite variance (post correction)', las=1, pch=16, log='xy', cex.lab=2)
# the residuals of run order on each metabolite, preserves the variation within each metabolite (unless there was runorder effects to remove).  It also shrinks the variance within each sample, because the metabolites are now all represented by residuals (mean=0).  I wonder if I were to remove the slope, but not the intercept from the run order lm?  Doesn't seem necessary, and does not acknowldge the differential sensitivity acorss metabolites of the LCMS

d[ ,mzs] <- tmp[ ,mzs] # accept the run order residuals

par(mfrow=c(3,4))
par(mar=c(2,2,2,2))
for (i in 1:12) {
  plot(d[ ,mois[i]] ~ d$Sample.order, pch=16, cex=0.5, xlab='run order', ylab='', las=1, main=mois[i])
  abline(lm(d[ ,mois[i]] ~ d$Sample.order)) }

for(i in 1:7) {
  plot(prcomp(d[ ,mzs], scale=T)$x[ ,i] ~ d$Sample.order, col=2, pch=16, main=paste('PC', i))} 

save(d, mzs, mz.info, vials, file='normalized LCMS data')
######################################




rm(list=ls())
######################################
load('normalized LCMS data')
# merge LCMS data with biological variables
d[1:6,1:10]
head(vials)

# age at each separation
sepDay1 <- 29 # 1st separation, 6/24, adjusted (3 added) day ~ 28.9)
sepDay2 <- 43 # 2nd separarion 7/8, adjusted day ~ 42.9
sepDays <- c(sepDay1, sepDay2) # age at each separation
d <- merge(vials, d)
table(d$Climbing)

d$Climbing[d$Climbing %in% 3:4] <- 'middle'
d$Climbing[d$Climbing == 1] <- 'top'
d$Climbing[d$Climbing == 2] <- 'mid-top'
d$Climbing[d$Climbing == 5] <- 'mid-bottom'
d$Climbing[d$Climbing == 6] <- 'bottom'

d$Climbing <- as.factor(d$Climbing)
levels(d$Climbing)
d$Climbing <- factor(d$Climbing, levels=c('top', 'mid-top', 'middle', 'mid-bottom', 'bottom'))
levels(d$Climbing)
table(d$Climbing)

d[1:6,1:10]
str(d)
d$Separation <- as.factor(letters[as.numeric(as.factor(paste(d$SepDate, d$Separation)))]) # make a unique name for each batch of flies run in the separation device
d$Separation
table(d$Separation, d$SepDate) # separation should completely segregate by sepdate
d$SepDate <- as.factor(d$SepDate)


save(d, mzs, mz.info, sepDays, file='LCMS_and_biological_data')
###################################################################################
rm(list=ls())
load('LCMS_and_biological_data')

levels(d$Climbing)
levels(d$SepDate) <- c('week 4', 'week 6')

######################################################
# analyse PCA
######################################################
pca <- prcomp(d[ ,mzs], scale=T)

eigs <- pca$sdev^2
tw <- AssocTests::tw(eigs, eigenL=length(eigs), criticalpoint =  2.0234)
sigEigs <- tw$SigntEigenL

propexp <- eigs/sum(eigs)
sum(propexp[1:4])

mod1 <- clmm(Climbing ~ SepDate + (1|Separation), d)
mod2 <- clmm(Climbing ~ SepDate + pca$x[ ,'PC3'] + (1|Separation), d)
mod3 <- clmm(Climbing ~ SepDate * pca$x[ ,'PC3'] + (1|Separation), d)
anova(mod1, mod2, mod3) # a main effect (ie. independent of age) of PC3 on the climbing rank is 'off the charts' P=5x10-6

pairs(cbind(d$Climbing, pca$x[ ,1:sigEigs]), col=as.numeric(d$Climbing), pch=19)
pairs(cbind(d$Climbing, pca$x[ ,1:sigEigs]), col=as.numeric(d$SepDate), pch=19)

pc <- as.data.frame(pca$x[ ,1:sigEigs])

# look for outliers
MD <- mahalanobis(pc, colMeans(pc), cov(pc))

hist(MD)

ggplot(d, aes(Climbing, MD, fill=Climbing, color=Climbing))+
  geom_boxplot(alpha=0.5)+
  geom_point()+
  labs(y=expression('Mahalanobis distance '*PC[1-4]))+
  theme_bw(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate)

which(MD>15) # sample 70 looks exceptional
pairs(pc[ ,1:4], col=ifelse(MD > 15, 'red', 'grey'), pch=19) # seems exceptional primarily among PC3, but there it is quite different, I will remove it
d <- d[-(which(MD>15)), ] # remove sample 70


save(d, mzs, mz.info, sepDays, file='LCMS_and_biological_data')
####################################################################################
rm(list=ls())

load('LCMS_and_biological_data')

pca <- prcomp(d[ ,mzs], scale=T) # recompute PCs
eigs <- pca$sdev^2
tw <- AssocTests::tw(eigs, eigenL=length(eigs), criticalpoint =  2.0234)
sigEigs <- tw$SigntEigenL

pc <- data.frame(pca$x[ ,1:sigEigs])

sumsq <- as.data.frame(cbind(
  anova(lm(pc$PC1 ~ d$SepDate * d$Climbing + d$Separation))$`Sum Sq`,
  anova(lm(pc$PC2 ~ d$SepDate * d$Climbing + d$Separation))$`Sum Sq`,
  anova(lm(pc$PC3 ~ d$SepDate * d$Climbing + d$Separation))$`Sum Sq`,
  anova(lm(pc$PC4 ~ d$SepDate * d$Climbing + d$Separation))$`Sum Sq` ))

rownames(sumsq) <- rownames(anova(lm(pc$PC1 ~ d$SepDate * d$Climbing + d$Separation)))
colnames(sumsq) <- paste0('PC', 1:ncol(sumsq))
colSums(sumsq)

sumsq
sumsq$variable <- c('age', 'climbing', 'replicate', 'age x climbing', 'residual')
sumsq <- gather(sumsq, key='PC', value='sumsq', -variable)
head(sumsq)
sumsq$PC <-  as.numeric(as.factor(sumsq$PC))

#############################
# a little scratch analysis:
anova(lm(pc$PC2 ~ d$SepDate * d$Climbing + d$Separation))
anova(lm(pc$PC3 ~ d$SepDate * d$Climbing + d$Separation))
summary(clm(Climbing ~ SepDate * pca$x[ ,'PC3'], data=d) ) # ordinal regression
anova(lm(pc$PC4 ~ d$SepDate * d$Climbing + d$Separation))
summary(clm(Climbing ~ SepDate * pca$x[ ,'PC4'], data=d) ) # ordinal regression

MD <- mahalanobis(pc, colMeans(pc), cov(pc)) # how about mehalonobis distance
anova(lm(MD ~ d$SepDate * d$Climbing + d$Separation))
summary(clm(Climbing ~ SepDate * MD, data=d) )

plot(MD ~ d$Climbing)
ggplot(d, aes(Climbing, MD, color=Climbing, group=Climbing, fill=Climbing))+
  geom_boxplot(alpha=0.5)+
  geom_point()+
  labs(y=expression('Mahalanobis distance '*PC[1-4]))+
  theme_bw(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate)

##############################

ggplot(sumsq, aes(fill=variable, y=sumsq, x=PC)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 14)+
  xlab('')+
  ylab('ANOVA sum of squares')+
  scale_fill_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "grey88", 'grey60'))

tmp <- d
tmp$SepDate
tmp$Age <- factor(tmp$SepDate)

tmp$PC1 <- pca$x[ ,'PC1']
tmp$PC2 <- pca$x[ ,'PC2']
tmp$PC3 <- pca$x[ ,'PC3']
tmp$PC4 <- pca$x[ ,'PC4']

tmp2 <- plyr::ddply(tmp, c("Age", "Climbing"), summarise, meanPC2=mean(PC2), meanPC3=mean(PC3))
head(tmp2)
tmp2 <- merge(tmp, tmp2)

scales::hue_pal()(12) # shows pallet used by ggplot with (n) colors

propexp <- eigs/sum(eigs)
sum(propexp[1:4]) # proportion of variance explained by significant PCs

sig <- 1:ncol(pca$x) %in% 1:sigEigs


Var <- data.frame('PropVarExp' = eigs/sum(eigs), 'PC'=colnames(pca$x), 'significant'=sig) 
Var$PC <- factor(Var$PC, levels=Var$PC)

head(Var)
levels(Var$PC)  <- 1:71

p1 <- ggplot(Var[1:10, ], aes(PC, PropVarExp, fill = significant))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c('grey88', "#00C08B"))+
  theme_bw(base_size = 16)+
  ylab('variance explained (%)')+
  theme(legend.position="none", axis.text.x = element_text(size = 8), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p2 <- ggplot(sumsq, aes(fill=variable, y=sumsq, x=PC)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 16)+
  ylab('ANOVA sum of squares')+
  scale_fill_manual(values=c("#F8766D", "orange", "#00BF7D", "grey88", 'grey60'))

ggarrange(p1, p2, widths = c(0.65, 1))


p3 <- ggplot(tmp2, aes(Age, meanPC2, color=Climbing, group=Climbing))+
  geom_line(alpha=0.8)+
  geom_point(aes(Age, PC2), position = position_jitter(w = 0.02, h = 0))+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  ylab('PC2')

p34 <- ggplot(tmp, aes(PC3, PC4, color=Climbing, group=Climbing, fill=Climbing))+
  geom_point(alpha=1)+
  stat_ellipse(geom = "polygon", alpha = 0.25, level=0.5)+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme(axis.text.x=element_blank())

ggarrange(p3, p34, common.legend = T, legend='right', widths = c(0.7, 1.2), nrow=1)

top_row <- ggarrange(p1, p2, widths = c(0.65, 1))
bot_row <- ggarrange(p3, p34, common.legend = T, legend='right', widths = c(0.7, 1.2), nrow=1)
ggarrange(top_row, bot_row, ncol=1)

# grab the legend from this plot to clarify the one on the bottom row (above)
ggplot(tmp, aes(PC3, PC4, color=Climbing, group=Climbing, fill=Climbing))+
  geom_point(alpha=1)+
  stat_ellipse(geom = "polygon", alpha = 0.25, level=0.5)+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme(axis.text.x=element_blank())




#####################################################################################
# univariate analysis
# Climbing is a complex variable. The climbing groups are not linearly related.  Try ordinal regression
str(d)

# ordinal regression, logisic model looking for metabolites that influence the odds of being a high or low ranked climber:
AICs <- matrix(nr=length(mzs), nc=3)
Ps <- matrix(nr=length(mzs), nc=2)

for(i in 1:length(mzs)) {
nullmod <- clmm(Climbing ~ SepDate + (1|Separation), data=d) 
mod1 <- clmm(Climbing ~ SepDate + d[ ,mzs[i]] + (1|Separation), data=d) 
mod2 <- clmm(Climbing ~ SepDate * d[ ,mzs[i]] + (1|Separation), data=d)  
a <- anova(nullmod, mod1, mod2)

AICs[i, ] <- a$AIC
Ps[i, ] <- a$`Pr(>Chisq)`[2:3]
}

colnames(Ps) <- c('age + mz', 'age x mz')
rownames(Ps) <- mzs

colnames(AICs) <- c('null', 'age + mz', 'age x mz')
rownames(AICs) <- mzs
AICs <- as.data.frame(AICs)

head(AICs[order(AICs$`age + mz`), ])
head(AICs[order(AICs$`age x mz`), ])

minAICs <- apply(AICs, 1, which.min) 

mzs[minAICs==3] # metabolites whos most explanatory model involves mz x age interaction on odds of climbing
mzs[minAICs==2] # metabolites whos most explanatory model involves a main effect of mz on odds of climbing

FDRs <- apply(Ps, 2, function(x) p.adjust(x, 'fdr'))
apply(FDRs, 2, function(x) table(x <=0.05)) # the non-interaction model was more explanatory than the interaction model for all 12 metabolties.  No evidence for mz ~ climbign x age, the effect of mz~climbing was not different acorss age.

hits <- mzs[rowSums(FDRs <=0.05)>0]
hits


save(hits, FDRs, file='ordinal_main_effect_analysis')



rm(list=ls())
#################################################################
load('LCMS_and_biological_data')
load('ordinal_main_effect_analysis')
# for plotting:
# get residuals of mz ~ separation (this will portray the mz values fit in the ordinal modeling above)
resids <- matrix(nr=nrow(d), nc=length(hits))

for(i in 1:length(hits)){
  resids[, i] <- lm(d[ ,hits[i]] ~ Separation, d)$residuals
}

colnames(resids) <- hits
tmp <- data.frame(d[ ,c('SepDate', 'Climbing')], resids) # this adds leading Xs to colnames with leading numbers
head(tmp)

# order hits and fix labels
hit.fdrs <- FDRs[hits,1]
hit.fdrs
names(hit.fdrs[order(hit.fdrs)])
names(hit.fdrs) <- colnames(tmp)[-c(1:2)]


# make a long df of the hits
long <- pivot_longer(tmp, cols=all_of(colnames(tmp)[-c(1:2)]), names_to = 'mz')
levels(long$SepDate) <- c('week 4', 'week 6')

table(long$mz)
# make metabolite names that are better for plot labels:
long$mz <- tolower(long$mz)

long$mz [long$mz=='alpha.ketoglutaric.acid'] <- 'a-ketogluterate'
long$mz [long$mz=='xanthurenic.acid'] <- 'xanthurenate'
long$mz [long$mz=='glycerol.3.p'] <- 'glycerol-3-P'
long$mz [long$mz=='linoleic.acid'] <- 'linolate'
long$mz [long$mz=='cgmp'] <- 'cGMP'
long$mz <- as.factor(long$mz)

levels(long$mz)
hits <- levels(long$mz)

p1 <- ggplot(subset(long, mz %in% levels(long$mz)[1:3]), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_point(size=1)+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  theme_bw()+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme(text = element_text(size=12), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')

p2 <- ggplot(subset(long, mz %in% levels(long$mz)[4:6]), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_point(size=1)+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme_bw()+
  theme(text = element_text(size=12), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')

p3 <- ggplot(subset(long, mz %in% levels(long$mz)[7:9]), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_point(size=1)+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme_bw()+
  theme(text = element_text(size=12), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')

ggarrange(p1, p2, p3, common.legend = T, legend='right', nrow=1)



ggplot(subset(long, mz %in% levels(long$mz)), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_point(size=1)+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme_bw(base_size = 14)+
  theme(text = element_text(size=9), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')


## make separate panels for mzs that are 'up' or 'down' among the bottom climbers:
hits
updown <- factor(c('d', 'u', 'u', 'd', 'u', 'u', 'd', 'd', 'u'))
names(updown) <- hits
table(updown)
long$updown <- updown[long$mz]
table(long$updown)

dplot <- ggplot(subset(long, mz %in% names(updown)[updown=='d']), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  geom_point(size=1)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme_bw(base_size = 14)+
  ylab('')+
  theme(text = element_text(size=14), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')

uplot <- ggplot(subset(long, mz %in% names(updown)[updown=='u']), aes(x=Climbing, y=value, group=SepDate, color=Climbing))+
  geom_smooth(method='loess', fill = "grey90", linewidth=0.5, color='grey50')+
  geom_point(size=1)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  theme_bw(base_size = 14)+
  ylab('')+
  theme(text = element_text(size=14), axis.text.x=element_blank()) +
  facet_grid(~mz ~ SepDate, scales='free')

ggarrange(dplot, uplot, common.legend = T, legend='right', nrow=1)




###########################################################################################################
# Enrichment analysis: FELLA
###########################################################################################################
suppressPackageStartupMessages({
  library(igraph)
  library(FELLA)
  library(org.Dm.eg.db)
  library(KEGGREST)
  library(magrittr)
  library(biomaRt)
  library(rvest)
  library(XML)
})

set.seed(1)

# get KEGG id for metabolites
head(mz.info)
hits

kegg <- mz.info$`KEGG ID`
hits.kegg <- mz.info$`KEGG ID`[mz.info$COMPOUND %in% hits]

# get compounds from the analysis
all_compounds <- mz.info$`KEGG ID`
all_compounds <- all_compounds[!is.na(all_compounds)]
all_compounds <- sapply(strsplit(all_compounds, "/"), '[[', 1)

hits.kegg <- mz.info$`KEGG ID`[mz.info$COMPOUND %in% hits]
hits.kegg <- hits.kegg[!is.na(hits.kegg)]
hits.kegg <- sapply(strsplit(hits.kegg, "/"), '[[', 1)

all_compounds
hits.kegg
getwd()
tmpdir <- "FELLA_database"

#############################################################################################################
# [IF YOU ALREADY RAN THE LINES HERE, SKIP TO THE LINE: fella.data <- loadKEGGdata ]
graph <- buildGraphFromKEGGREST(organism = "dme", filter.path = c("01100", "01200", "01210", "01212", "01230"))  # Build network and filter out the large pathways:
unlink(tmpdir, recursive = TRUE) # this line allows you to rewrite the: buildataFromGraph()

buildDataFromGraph(graph, databaseDir = tmpdir, internalDir = FALSE, matrices = 'diffusion', normality = "diffusion", niter = 1000) # using listMethods() here leads to three matrices being made: 'diffusion', ' hypergemetric' and 'pagerank'.  It takes longer than building with a single method.  Note that 'normality' can either be 'diffusion' or 'pagerank', and this will determine which of those two analyses can run based on these data.  I ran normality as 'diffusion' and this caused an error saying that Z-scores would not be available for the 'pagerank' method.  It did however make a 'pagerank' matrix, so I suspect that I can still do non-parametric analysis by pagerank, rather than a parametric analysis that might depend on the outcome of this line.
#############################################################################################################

fella.data <- loadKEGGdata(databaseDir = tmpdir, internalDir = FALSE, loadMatrix = 'diffusion') # using listMethods() here loads all 3 matrices, which in turn allows all 3 analysis methods
getInfo(fella.data) # prints the version KEGG used (Sergio says it always uses the most recent release)
fella.data

## enrichment analysis
enrichment_analysis <- enrich(compounds = hits.kegg, compoundsBackground = all_compounds, data = fella.data, method = 'diffusion', approx = "simulation", niter=10000)
enrichment_table <- generateResultsTable(object=enrichment_analysis, data=fella.data)
enrichment_table$FDR <- p.adjust(enrichment_table$p.score, 'BH')
head(enrichment_table[order(enrichment_table$p.score), ]) 
# find and replace commas with underscores to enable r to read the results later
enrichment_table$KEGG.name <- gsub(',', '_', enrichment_table$KEGG.name)

enrichment_table[enrichment_table$Entry.type %in% c('pathway', 'module'), ] # a look at pathways and modules only

write.table(enrichment_table, 'FELLA_results/FELLA_table.csv', sep=',', row.names = F, quote=F)


################################################################################################
dir('FELLA_results')

# [NOTE]: open enrichment tables in excel, find and replace ',' with '_', then save as .csv, then you can open with:
fella <- read.csv('FELLA_results/FELLA_table.csv') 
head(fella)
table(fella$Entry.type)

# limit to sig figs:
fella[ ,4:5] <- apply(fella[ ,4:5], 2, function(x) signif(x, 3))
# simplify kegg names
fella$KEGG.name <- sapply(str_split(fella$KEGG.name, "\\ -"), '[[', 1) 
fella$KEGG.name <- sapply(strsplit(fella$KEGG.name, split='_'), '[[', 1) 

DT::datatable(fella) # save an HTML table - may not properly identify delimiters


# create static table for publication
library(gt)
library(gapminder)

fella.table.fig <- fella %>% arrange(p.score) %>% filter(Entry.type %in% c('pathway', 'module')) %>% gt(groupname_col = "Entry.type", rowname_col = "KEGG.name") %>% tab_header(title = "Table 3. Pathway Enrichment") %>% cols_align(align = "center") %>% row_group_order(groups = c('pathway')) # [NOTE]: add 'module' into the last argument to list any modules identified.

fella.table.fig

dir('FELLA_results')
save(file='FELLA_results/fella.data', fella, fella.data, all_compounds, hits.kegg, hits, mz.info, enrichment_analysis)

rm(list=ls())
################################################################################################



