#August 8, 2013
#Function to create life table calculations
#Scott Pletcher spletch@umich.edu
suppressPackageStartupMessages({
  library(splitstackshape)
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  
})

rm(list=ls())

setwd("/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/Claudia Sun")
dir()

# load Scott's Life Table Function
source('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/R.functions/LT.Pletcher.R')

# [Summary of the study]: From May to August 2021, Claudia Sun ran an experiment with Canton-S mated females.  5/29/2021 was the first recorded census. Flies were likely 2-4 days old at that point (*see 3 days (108h) added to ageH, below*) . Flies were raised in density controlled bottles and were on Promilsow lab 'standard food' the whole time.

# there were 132 vials set up in total, each with 25 females at t=0
# on two different days, 6/24/2021 (week~4, day=31) and 7/8/2021 (week~6, day=45), some vials were removed from the demography, pooled, and Claudia ran 'separations', where batches of pooled flies (4 batches on 6/24, and 7 batches on 7/8, see below) were put into a polystyrene tube (dimensions= __ x__), banged to the mid-bottom, allowed to climb until ~1/2 the flies made it across approximately the middle of the tube, then a paper separator was slid into the middle of the tube and the flies on the mid-top and on the mid-bottom were poured into new vials.  This process was repeated 2 more times for the flies initially reaching the mid-top, and 2 more times for those on the mid-bottom.  This left 6 fractions, or levels of Climbing behavior.  These Climbing groups were then anesthetized on light CO2 (Ben did this), split into ~20-30 flies per vial, and these new vials were put back into the demography.  At each separation there were fewer vials returned to the demography than there were vials removed from the demography (see **Vial Assignments** below).   

# read descriptions of the vials:
vials <- read.table("Claudias data and code/Vial-Specific info.csv", sep= ",", head = T)
head(vials)

# **Vial Assignments**:
# vials 1-40 were pooled for the 1st separation, which gave vials 1-32 after the separation.
# vials 41-124 were pooled for the 2nd separation, which gave vials 57-105 after the separation.
# vials 125-132 were never removed from the demography and represent unmolested survivorship

table(vials$SepDate) # the number of vials that were generated after each separation day was complete
table(vials$SepDate, vials$Separation) # the batches of flies that were separated together have strange named (Ben's fault).

dat <- read.table("Claudias data and code/climbers_lifespan_A.csv", sep= ",", head = T)
head(dat)
added.age <- 3*24
dat$AgeH <- dat$AgeH + added.age


# a few observations from Claudia's notebook (pg 14):
# 7/6/21 census: all of the flies in vial 88 were flipped into a new vial of food that was wet, all of these flies were killed
dat[grepl('7/', dat$CensusTime) & dat$Chamber==88, ] # on 7/8/21 (the next census interval), this vial was replaced with flies from the 2nd separation. There is no need to do any correction to the data as they are.

# 7/8/21 census: many of the flies in vials 70-90 died likely because the vials they were flipped into on 7/6/21 had water in them.
dat[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90, ] # these IntDeaths should instead be made Censored, and replace with 0
dat$Censored[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90] <- dat$Censored[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90] + dat$IntDeaths[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90]
dat$IntDeaths[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90] = 0
dat[grepl('7/8/21', dat$CensusTime) & dat$Chamber%in% 70:90, ] 

# *NOTE: the previous two points were observed for flies prior to the separation.  The next observation was made on flies that had been separated:*

# 7/10/21 census: vial 8 had a lot of water, likely causing drowning
dat[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8, ] # IntDeaths flies for this vial in this census should be added to Censored and replaced with 0
dat$Censored[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8] <- dat$Censored[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8] + dat$IntDeaths[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8]
dat$IntDeaths[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8] = 0
dat[grepl('7/10/21', dat$CensusTime) & dat$Chamber==8, ] 


x <- merge(vials, dat, by.y= "Chamber", by.x= "vial")
head(x)
x$Climbing <- as.factor(x$Climbing)
x$day <- x$AgeH/24
x$SepDate[x$vial %in% 125:132] <- 'demography.control'
x$SepDate[x$vial %in% 33:40] <- 'removed_6/24/2021' 
x$SepDate[x$vial %in% c(41:56, 106:124)] <- 'removed_7/8/2021'
table(x$SepDate)
head(x)

# 1st separation, 6/24, adjusted (4.5 added) day = 30.44292)
# 2nd separarion 7/8, adjusted day = 44.359167

# since there is no way to know (beyond estimation) what was the population size in the vials used in each separation prior to the separation, I cannot know the at-risk population in them.  So, the data from vials other than the 'demography.control' should not be used for mortality analysis.
table(x$SepDate == '6/24/2021' & x$day > 30.5) # N observations from sep 1
table(x$SepDate == '7/8/2021' & x$day > 44.4)  # N observations from sep 2
table(x$SepDate=='demography.control') # N observations from the control population
table((x$SepDate == '6/24/2021' & x$day > 30.5)|(x$SepDate == '7/8/2021' & x$day > 44.4)|(x$SepDate=='demography.control')) # total flies observed

presep <- x[x$SepDate=='demography.control', ] # collect all data for flies not subject to separation
presep$population <- 'demography.control'
presep$Climbing = NA
presep$Separation = NA
presep$SepDate = NA

postsep <- x[((x$SepDate == '6/24/2021' & x$day > 30.5) | (x$SepDate == '7/8/2021' & x$day > 44.4)), ] # collect all flies that were in sep1 or sep2
table(postsep$SepDate)
postsep$population <- ifelse(postsep$SepDate == '6/24/2021', 'week4', 'week6')

dat <- rbind(presep, postsep)
head(dat)
table(dat$SepDate)
dat$SepDate <- as.factor(dat$SepDate)
dat$population <- as.factor(dat$population)

table(dat$Climbing)
dat$Climbing[dat$Climbing %in% c('3', '4')] <- 'middle' 
dat$Climbing[dat$Climbing == '5'] <- 'mid-bottom'
dat$Climbing[dat$Climbing == '6'] <- 'bottom'
dat$Climbing[dat$Climbing == '1'] <- 'top'
dat$Climbing[dat$Climbing == '2'] <- 'mid-top'
dat$Climbing <- as.factor(dat$Climbing)
levels(dat$Climbing)
dat$Climbing <- factor(dat$Climbing, levels=c('top', 'mid-top', 'middle', 'mid-bottom', 'bottom'))
levels(dat$Climbing)

table(dat$Climbing, dat$population) 
table(dat$SepDate , dat$population)
table(dat$Separation , dat$population)

sepRep <- as.factor(paste(dat$SepDate, dat$Separation))
table(sepRep)
levels(sepRep) <- c(letters[1:11], NA)
dat$sepRep <- sepRep

# mortality analysis can be done using Scott Pletcher's life table-generating function: LT(), which needs a 3-column numeric matrix:
# create a numeric matrix of Three columns (time, number of events, event).
## event = 0 is censored (or alive) and event = 1 is dead.

n.events <- c(dat$Censored, dat$IntDeaths)
event <- c(rep(0, length(dat$Censored)), rep(1, length(dat$IntDeaths)) )
time <- rep(dat$day, 2)
mat <- cbind(time, n.events, event)
dim(mat)

# experimental variables
population <- rep(dat$population, 2) 
Climbing <- rep(dat$Climbing, 2) 
sepRep <- rep(dat$sepRep, 2) 
sepWeek <- rep(dat$SepDate, 2)

cond <- data.frame(population=population, Climbing=Climbing, sepWeek=sepWeek, sepRep=sepRep)
head(cond)

LTdat <- data.frame(cond, mat)
head(LTdat)
str(LTdat)


min(LTdat$time[LTdat$population == 'week6'])
min(LTdat$time[LTdat$population == 'week4'])

# loop thorough conditions, making life tables for each 
mort <- NULL

pops <- unique(population)

for(k in 2:3) {
  temp <- LTdat[LTdat$population== pops[k], c(2,5:7)]
  head(temp)
  trts <- levels(droplevels(temp$Climbing))
  
  for(i in 1:length(trts)) {
    x <- as.matrix(temp[temp$Climbing == as.character(trts[i]), 2:4]) 
    ltx <- LT(x)
    ltx$Climbing <- trts[i]
    ltx$population <- pops[k]
    mort <- rbind(mort, ltx)
  }}

head(mort)

# get mortality for the base population
x <- as.matrix(LTdat[LTdat$population == 'demography.control', 5:7]) 
ltx <- LT(x)
ltx$Climbing <- NA
ltx$population <- 'demography.control'
mort <- rbind(mort, ltx)

head(mort)

mort$Climbing[is.na(mort$Climbing)] <- 'control'

str(mort)
mort$Climbing <- as.factor(mort$Climbing)
levels(mort$Climbing)
mort$Climbing <- factor(mort$Climbing, levels=c('control', levels(dat$Climbing)))
table(mort$Climbing, mort$population)

mort$population <- as.factor(mort$population)

table(is.finite(mort$lnux2))
mort <- mort[is.finite(mort$lnux2), ]

head(mort)


######################################################################  
# probability of survival up to the first and second separations:
prob.surv.till.1st.sep <- 1-sum(mort$qx[mort$population == 'demography.control' & mort$IntEnd < 29])
flies.per.1st.sep <- 8*25
flies.per.1st.sep <- flies.per.1st.sep * prob.surv.till.1st.sep

prob.surv.till.2nd.sep <- 1-sum(mort$qx[mort$population == 'demography.control' & mort$IntEnd < 43])
flies.per.2nd.sep <- 7*25
flies.per.2nd.sep <- flies.per.2nd.sep * prob.surv.till.2nd.sep
######################################################################  
save(dat, mort, file='lifetableResults')


rm(list=ls())
###################################################################### 
# control population plot, to show context of the flies removed for fractionation
load('processed.survival.data.Claudias.expt')
par(mfrow=c(1,1))
s <- survfit(Surv(day, event) ~1, d[d$SepDate == 'demography.control', ])
plot(s, conf.int=F, lwd=2, las=1, xlab='time (day)', ylab='survivorship', main='demography control, n=193 flies, 8 vials')
abline(v=29)
abline(v=43)

par(mfrow=c(2,1))
plot(s, conf.int=F, lwd=2, las=1, xlab='time (day)', ylab='survivorship', col='black')
hist(d$day[d$SepDate == 'demography.control'], border=0, 40, col='grey30', main='', xlab='time of death (day)', ylab='n flies', las=1)
###################################################################### 


rm(list=ls())

######################################################################  
load('lifetableResults')
load('processed.survival.data.Claudias.expt') # loads d, which is the survival data from Claudia's demography, and 'sepDays', which are the days the segregations were done
load('gompertz.alphas') # load 'gfits', which are the alpha and beta of Gompertz fit with ML (in flexsurv)

head(d)


mort$Climbing
mort <- subset(mort, !Climbing %in% 'control')


alphas
mort$population

p3 <- ggplot(subset(mort, population=='week4'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point(pch=20) +
  theme_classic(base_size = 14)+
  scale_color_brewer(palette = "PuOr")+
  ggtitle('week 4 fractionation')+
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week4')
  )

p4 <- ggplot(subset(mort, population=='week6'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point(pch=20) +
  theme_classic(base_size = 14)+
  scale_color_brewer(palette = "PuOr")+
  ggtitle('week 6 fractionation')+
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week6'))
ggarrange(p3, p4, common.legend = T, legend='right')


ggplot(subset(mort, population=='week4'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point() +
  theme_bw(base_size = 14)+
  scale_color_brewer(palette = "PuOr")+
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week4')
  )


# some alternative plotting strategies:
middle <- subset(mort, Climbing=='middle')

ggplot(subset(mort, population=='week4' & Climbing != 'middle'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point() +
  theme_bw(base_size = 14)+
  facet_wrap(~Climbing)+
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week4' & Climbing != 'middle')
  )

ggplot(data=subset(mort, population=='week4' & Climbing != 'middle'), aes(y=lnux2, x=IntMid)) + 
  geom_point(data=select(middle,-Climbing), colour="grey80") +
  geom_point(size=2, pch=19, aes(color=Climbing)) +
  facet_wrap(~Climbing) +
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week4' & Climbing != 'middle'))+
  theme_bw(base_size = 14)

ggplot(data=subset(mort, population=='week4' & Climbing != 'middle'), aes(y=lnux2, x=IntMid)) + 
  geom_point(data=select(middle,-Climbing), colour="grey80") +
  geom_point(size=2, pch=19, aes(color=Climbing)) +
  facet_wrap(~Climbing) +
  geom_abline(mapping = aes(intercept = alpha, slope = beta, col = Climbing), data = subset(alphas, SepDate=='week4' & Climbing != 'middle'))+
  geom_abline(mapping = aes(intercept = alpha, slope = beta, color='grey80'), data = subset(alphas, SepDate=='week4' & Climbing == 'middle'))+
  theme_bw(base_size = 14)









