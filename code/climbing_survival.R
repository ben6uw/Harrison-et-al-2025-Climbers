
library(splitstackshape)
library(survival)
library(ggplot2)
library(survminer)
library(flexsurv)


rm(list=ls())

setwd("/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/Claudia Sun")
dir()

# [Summary of the study]: From May to August 2021, Claudia Sun ran an experiment with Canton-S mated females.  5/29/2021 was the first recorded census. Flies were likely 2-4 days old at that point (*see 3 days (72h) added to ageH, below*) . Flies were raised in density controlled bottles and were on Promilsow lab 'standard food' the whole time.

# there were 132 vials set up in total, each with 25 females at t=0
# on two different days, 6/24/2021 (week~4, day=31) and 7/8/2021 (week~6, day=45), some vials were removed from the demography, pooled, and Claudia ran 'separations', where batches of pooled flies (4 batches on 6/24, and 7 batches on 7/8, see below) were put into a polystyrene tube (dimensions= __ x__), banged to the mid-bottom, allowed to climb until ~1/2 the flies made it across approximately the middle of the tube, then a paper separator was slid into the middle of the tube and the flies on the mid-top and on the mid-bottom were poured into new vials.  This process was repeated 2 more times for the flies initially reaching the mid-top, and 2 more times for those on the mid-bottom.  This left 6 fractions, or levels of climbing behavior.  These climbing groups were then anesthetized on light CO2 (Ben did this), split into ~20-30 flies per vial, and these new vials were put back into the demography.  At each separation there were fewer vials returned to the demography than there were vials removed from the demography (see **Vial Assignments** below).   

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


sum(dat$IntDeaths) # total observed deaths
sum(dat$Censored) # total observed censored flies (eg. flew away, crushes or were stuck)

d <- merge(vials, dat, by.y= "Chamber", by.x= "vial")
head(d)

junk <- c('Collective', 'UniqueName', 'Transfer', 'Flag1', 'Flag2', 'N', 'Carried', 'Deaths')
d <- d[d$AgeH != 1, ] # remove first 'dummy' census (an artifact of Dlife)
d <- d[ , !colnames(d) %in% junk]
head(d)

d$Climbing[d$Climbing %in% 3:4] <- 'middle'
d$Climbing[d$Climbing == 1] <- 'top'
d$Climbing[d$Climbing == 2] <- 'mid-top'
d$Climbing[d$Climbing == 5] <- 'mid-bottom'
d$Climbing[d$Climbing == 6] <- 'bottom'

d$Climbing <- as.factor(d$Climbing)
levels(d$Climbing)
d$Climbing <- factor(d$Climbing, levels=c('top', 'mid-top', 'middle', 'mid-bottom', 'bottom'))
levels(d$Climbing)

length(unique(d$CensusTime)) # 47 (real) censuses
d$day <- d$AgeH/24
head(d)

# **Vial Assignments**:
# vials 1-40 were pooled for the 1st separation, which gave vials 1-32 after the separation.
# vials 41-124 were pooled for the 2nd separation, which gave vials 57-105 after the separation.
# vials 125-132 were never removed from the demography and represent unmolested survivorship

table(d$SepDate, d$vial)

colSums(is.na(d))# NA here means these vials were either the unremoved control vials (vials 125-132), or were vials that were not replaced after a separation had been done.
d$SepDate[d$vial %in% 125:132] <- 'demography.control'
d <- d[!d$vial %in% 33:40, ] # remove preseparation data
d <- d[!d$vial %in% c(41:56, 106:124), ]  # remove preseparation data
table(d$SepDate)
d$SepDate <- as.factor(d$SepDate)
levels(d$SepDate)
d$SepDate <- factor(d$SepDate, levels=c('demography.control', '6/24/2021', '7/8/2021'))
levels(d$SepDate)
levels(d$SepDate) <- c('demography.control', 'week4', 'week6')
# keep only data from vials from demography control, or those from the separations but only after the separations: the pre-separation data is not usable, particularly because the number of flies censused when the vials were pooled is not known (can only be estimated).

sepDay1 <- 29 # 1st separation, 6/24, adjusted (3 added) day ~ 28.9)
sepDay2 <- 43 # 2nd separarion 7/8, adjusted day ~ 42.9

# remove preseparation data
table((d$SepDate == 'week4' & d$day < 30) | (d$SepDate == 'week6' & d$day < 44))
d <- d[!((d$SepDate == 'week4' & d$day < 30) | (d$SepDate == 'week6' & d$day < 44)), ]

table(d$SepDate) # total remaining observations per population

# make 'events'

death <- expandRows(d, 'IntDeaths')
censored <- expandRows(d, "Censored")
death$event <- 1
censored$event <- 0

death <- subset( death, select = -Censored )
censored <- subset( censored, select = -IntDeaths)
d <- rbind(death, censored)
rm(death)
rm(censored)

table(d$event) # total recorded deaths and censored events
head(d)
str(d)

# distributions of deaths and censored flies over age
events <- ifelse(d$event == 1,'death', 'censored')
sepDays <- c(sepDay1, sepDay2)

save(d, sepDays, file='processed.survival.data.Claudias.expt')
##################################################################







rm(list=ls())
################################################################## 
load('processed.survival.data.Claudias.expt')
##################################################################

ggplot(d, aes(x = day, fill = SepDate)) + 
  geom_histogram(alpha = 0.7) +
  theme_classic() +
  labs(title="Event distribution in Claudia's experiment") +
  theme(plot.title = element_text(size=16), axis.text=element_text(size=14), axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  scale_fill_discrete(name = "")


par(mfrow=c(1,2))

table(d$SepDate) # total observations per population
colSums(table(d$vial, d$SepDate)>0) # number of vials per population
table(d$Separation, d$SepDate)
table(d$Climbing, d$SepDate) 
range(table(d$Climbing, d$SepDate)[table(d$Climbing, d$SepDate)>0]) # range of sample size for each climbing group

d$tmp <- 1
d <- d %>% group_by(vial) %>% mutate(n.flies = sum(tmp)) # add a column of N flies per vial after the fractionation 
d = subset(d, select = -c(tmp) )

x <- table(d$vial[d$SepDate !='demography.control'])
mean(x)
sd(x)
range(x)
hist(x, xlab='flies per vial, post separation', main='', border=0, col='grey40', 20)

# analyze flies per vial as a covariate on mortality
ggplot(subset(d, SepDate !='demography.control') , aes(AgeH/24, x=n.flies, color=Climbing))+
  geom_point()+
  theme_classic()+
  ylab('age at event (day)')+
  xlab('flies in vial after fractionation')+
  geom_smooth(method='lm')+
  ggh4x::facet_grid2(~SepDate~Climbing, scales = "free", independent = "all")
# the number of flies in a vial following fracitonation is not associated with remaining lifespan.


save(d, sepDays, file='dataForSurvivalModels')




rm(list=ls())
#############################################################################
load('dataForSurvivalModels')

s <- survfit(Surv(day, event) ~1, d[d$SepDate == 'demography.control', ])
plot(s, conf.int=F, lwd=2, las=1, xlab='time (day)', ylab='survivorship', main='demography control, n=193 flies, 8 vials')
abline(v=29)
abline(v=43)

# partition data into post-separation populations
sep1 <- d[d$day > 32 & d$SepDate == "week4", ]
sep1$SepDate <- droplevels(sep1$SepDate)
sep2 <- d[d$day > 45 & d$SepDate == "week6", ] 
sep2 $SepDate <- droplevels(sep2$SepDate)

library(RColorBrewer)
brewer.pal(n=5,"PuOr")
brewer.pal(n=5,"Spectral")

s1 <- survfit(Surv(day, event) ~ Climbing, sep1)
p1 <- ggsurvplot(s1, data=sep1, legend='right', legend.title="", xlim=c(sepDays[1], 100), size=0.5, censor.size=2, xlab='Age (day)') 
p1 <- p1$plot + scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))
p1

s2 <- survfit(Surv(day, event) ~ Climbing, sep2)
p2 <- ggsurvplot(s2, data=sep2, legend='right', legend.title="", xlim=c(sepDays[2], 100), size=0.5, censor.size=2, xlab='Age (day)') 
p2 <- p2$plot + scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))
p2

ggarrange(p1, p2, common.legend = T)


# get flexsurv fits of alphas and betas
tmpList <- list()

for(k in 2:3) {
  tmp <- d[d$SepDate ==levels(d$SepDate)[k], ]
  tmp$Climbing <- factor(tmp$Climbing, levels=c("middle", "top", "mid-top", "mid-bottom", "bottom"))
  gom <- flexsurvreg(Surv(day, event) ~ Climbing + shape(Climbing), data=tmp, dist = "gompertz")
  alphas <- as.data.frame(gom$res.t[3:6, 1:3] + gom$res.t[2,1])
  alphas 
  alphas$se <- gom$res.t[3:6, 4]
  gom$res.t[2, ] # middle climbers alpha 
  alphas <- rbind(gom$res.t[2, ], alphas) # alphas of all climbing groups
  alphas
  colnames(alphas) <- paste0('alpha_', colnames(alphas))
  
  betas <- as.data.frame(gom$res.t[7:10, 1:3] + gom$res.t[1,1])
  betas$se <-  gom$res.t[7:10, 4]
  gom$res.t[1, ] # middle climbers alpha 
  betas <- rbind(gom$res.t[1, ], betas) # alphas of all climbing groups
  colnames(betas) <- paste0('beta_', colnames(betas))
  betas
  goms <- data.frame(alphas, betas)
  goms$Climbing <- as.factor(levels(tmp$Climbing))
  goms$SepDate <- unique(tmp$SepDate) 
  tmpList[[k-1]] <- goms}

goms <- do.call(rbind, tmpList)


# compare models with and without beta:
tmp <- d[d$SepDate =='week4', ]
ab <- flexsurvreg(Surv(day, event) ~ Climbing + shape(Climbing), data=tmp, dist = "gompertz") # alpha and beta
a <- flexsurvreg(Surv(day, event) ~ Climbing, data=tmp, dist = "gompertz") # alpha only
null <- flexsurvreg(Surv(day, event) ~ 1, data=tmp, dist = "gompertz") # null (no fractions)
lmtest::lrtest(null, a, ab) # liklihood ratio test

tmp <- d[d$SepDate =='week6', ]
ab <- flexsurvreg(Surv(day, event) ~ Climbing + shape(Climbing), data=tmp, dist = "gompertz") # alpha and beta
a <- flexsurvreg(Surv(day, event) ~ Climbing, data=tmp, dist = "gompertz") # alpha only
null <- flexsurvreg(Surv(day, event) ~ 1, data=tmp, dist = "gompertz") # null (no fractions)
lmtest::lrtest(null, a, ab) # liklihood ratio test
# beta is singnificant after week 4, but not after week6
# intrinsic mortality (alpha) is singificant after fractionation at both ages


save(goms, file='gompertz.fits')
save(sep1, sep2, sepDays, file='sep.results')



gom <- flexsurvreg(Surv(day, event) ~ Climbing + shape(Climbing), data=tmp, dist = "gompertz")
gom



rm(list=ls())
##########################################################################################
load('processed.survival.data.Claudias.expt')
load('sep.results')
load('lifetableResults') # gets mortality estimates from Scotts life table function
load('gompertz.fits') # load 'gfits', which are the alpha and beta of Gompertz fit with ML (in flexsurv)

s1 <- survfit(Surv(day, event) ~ Climbing, sep1)
p1 <- ggsurvplot(s1, data=sep1, legend='right', legend.title="", xlim=c(sepDays[1], 100), size=0.5, censor.size=2, xlab='age (day)', ylab='survival probability') 
p1 <- p1$plot + scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))
p1

s2 <- survfit(Surv(day, event) ~ Climbing, sep2)
p2 <- ggsurvplot(s2, data=sep2, legend='right', legend.title="", xlim=c(sepDays[2], 100), size=0.5, censor.size=2, xlab='age (day)', ylab='survival probability') 
p2 <- p2$plot + scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))
p2

goms

p3 <- ggplot(subset(mort, population=='week4'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  xlab('age (day)')+
  geom_point(pch=20) +
  theme_classic(base_size = 16)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  geom_abline(mapping = aes(intercept = alpha_est, slope = beta_est, col = Climbing), data = subset(goms, SepDate=='week4')
              , alpha=0.8)

p4 <- ggplot(subset(mort, population=='week6'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  xlab('age (day)') +
  geom_point(pch=20) +
  theme_classic(base_size = 16)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  geom_abline(mapping = aes(intercept = alpha_est, slope = beta_est, col = Climbing), data = subset(goms, SepDate=='week6'), alpha=0.8)

ggarrange(p3, p4, common.legend = T, legend='right')

ggarrange(p1, p2, p3, p4, common.legend = T, legend='right')


head(goms)

alphaPlot <- ggplot(goms, aes(y=alpha_est, x=Climbing, group=SepDate, color=Climbing))+
  geom_point(pch=20, size=4) +
  theme_classic(base_size = 14)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  facet_wrap(~SepDate)+
geom_errorbar(aes(ymin=alpha_est-alpha_se, ymax=alpha_est+alpha_se), width=.2, position=position_dodge(.9)) 


betaPlot <- ggplot(goms, aes(y=beta_est, x=Climbing, group=SepDate, color=Climbing))+
  geom_point(pch=20, size=4) +
  theme_bw(base_size = 14)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  facet_wrap(~SepDate)+
  geom_errorbar(aes(ymin=beta_est-beta_se, ymax=beta_est+beta_se), width=.2, position=position_dodge(.9)) 


ggarrange(alphaPlot, betaPlot, common.legend = T, legend='right')



##################################################################
# plot mean lifespan +/- error per group
st1 <- survival:::survmean(s1, rmean=max(d$day)) # stats for climbing groups in sep 1
st2 <- survival:::survmean(s2, rmean=max(d$day)) # stats for climbing groups in sep 2

st1 <- as.data.frame(st1$matrix)
st2 <- as.data.frame(st2$matrix)

st1$SepDate='week4'
st2$SepDate='week6'

st1$Climbing <- as.factor(gsub('Climbing=', '', rownames(st1)))
st1$Climbing <- factor(st1$Climbing, levels=c('top', 'mid-top', 'middle', 'mid-bottom', 'bottom'))
st1$Climbing

st2$Climbing <- as.factor(gsub('Climbing=', '', rownames(st2)))
st2$Climbing <- factor(st2$Climbing, levels=c('top', 'mid-top', 'middle', 'mid-bottom', 'bottom'))
st2$Climbing

stall <- data.frame(rbind(st1, st2))

head(stall)

ggplot(stall, aes(y=rmean, x=Climbing, group=SepDate, color=Climbing))+
  geom_line(col='grey50')+
  geom_errorbar(aes(ymin=X0.95LCL, ymax=X0.95UCL), width=0.05, color='grey50')+
  geom_point(size=2)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate, scales='free')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  ylab('mean lifespan (days +/-95% CI)')

ggplot(stall, aes(y=rmean, x=Climbing, group=SepDate, color=Climbing))+
  geom_line(col='grey50')+
  geom_errorbar(aes(ymin=rmean-se.rmean., ymax=rmean+se.rmean.), width=0.05, color='grey50')+
  geom_point(size=2)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate, scales='free')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  ylab('mean lifespan (days +/-se)')



##################################################################
# combine with Harini's age-of-onset data:

original.stall <- stall

p1 <- ggplot(original.stall, aes(y=rmean, x=Climbing, group=SepDate, color=Climbing))+
  geom_line(col='grey50')+
  geom_errorbar(aes(ymin=rmean-se.rmean., ymax=rmean+se.rmean.), width=0.05, color='grey50')+
  geom_point(size=2)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate, scales='free')+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  ylab('mean lifespan (days +/-se)')

p1


load("~/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/Claudia Sun/Harini_Shankar/Climbing Lifespan/Harinis.data")
# mean lifespan +/- error per group from Harini's analysis
st1 <- survival:::survmean(survfit(Surv(age, event) ~ climbers, data=subset(dat, week==1)), rmean=max(dat$age)) 
st2.5 <- survival:::survmean(survfit(Surv(age, event) ~ climbers, data=subset(dat, week==2.5)), rmean=max(dat$age))  
st4 <- survival:::survmean(survfit(Surv(age, event) ~ climbers, data=subset(dat, week==4)), rmean=max(dat$age))

st1 <- as.data.frame(st1$matrix)
st2.5 <- as.data.frame(st2.5$matrix)
st4 <- as.data.frame(st4$matrix)

st1$SepDate='week1'
st2.5$SepDate='week2.5'
st4$SepDate='week4'

st1$Climbing <- as.factor(gsub('climbers=', '', rownames(st1)))
st2.5$Climbing <- as.factor(gsub('climbers=', '', rownames(st2.5)))
st4$Climbing <- as.factor(gsub('climbers=', '', rownames(st4)))

harini_stall <- data.frame(rbind(st1, st2.5, st4))

harini_stall
harini_stall$Climbing
harini_stall$Climbing <- factor(harini_stall$Climbing, levels=c('top', 'bot'))
levels(harini_stall$Climbing)[2] <- 'bottom'

p2 <- ggplot(harini_stall, aes(y=rmean, x=Climbing, group=SepDate, color=Climbing))+
  geom_line(col='grey50')+
  geom_errorbar(aes(ymin=rmean-se.rmean., ymax=rmean+se.rmean.), width=0.05, color='grey50')+
  geom_point(size=2)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~SepDate)+
  scale_color_manual(values=c("red","#5E3C99"))+
  ylab('mean lifespan (days +/-se)')


ggarrange(p1, p2, ncol=1, common.legend = T, legend='right')

p1


head(mort)
head(goms)


ggplot(subset(mort, population != 'demography.control'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point(pch=20) +
  theme_classic(base_size = 14)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  facet_wrap(~population)

ggplot(subset(mort, population != 'demography.control'), aes(y=lnux2, x=IntMid, group=Climbing, color=Climbing))+
  labs(x='Age (day)') +
  geom_point(pch=20) +
  geom_line(alpha=0.5)+
  theme_classic(base_size = 14)+
  scale_color_manual(values=c("red", "#FDB863", "grey", 'purple', "#5E3C99"))+
  facet_wrap(~population)



# does the data fit the PH assumption of cox?
# I doubt it, given the gomperts fits:
?coxph

tmp <- d[d$SepDate != 'demography.control', ] # fit dat from both ages
tmp$SepDate <- droplevels(tmp$SepDate)
tmp$Climbing <- factor(tmp$Climbing, levels=c('middle', 'bottom', 'mid-bottom', 'mid-top', 'top'))
head(tmp)
coxph(Surv(day, event) ~ Climbing * SepDate, data=tmp)
cox.zph(coxph(Surv(day, event) ~ Climbing * SepDate, data=tmp)) # violated

# try Aalen regression:
aaW4 <- aareg(Surv(day, event) ~ Climbing, tmp[tmp$SepDate=='week4', ])
aaW6 <- aareg(Surv(day, event) ~ Climbing, tmp[tmp$SepDate=='week6', ])
aaW4
aaW6

plot(aaW4, las=1)
abline(h=0, lty=2, col=2)

?aareg

# i tested the idea that left censorship was affecting the regression, by subtracting 30d from the week4 data, but the outcome of the aareg was identical

par(mar=c(4,4,2,2))
par(mfcol=c(5,2))
plot(aaW4, xlim=c(30, 85), las=1)
plot(aaW6, xlim=c(44, 85), las=1)


library(ggfortify)

w4 <- autoplot(aaW4, nrow=1, xlab='age (days)', ylab=expression(paste("hazard (", italic("H"), ")")))+
  geom_abline(intercept = 0, slope=0, linetype='dashed')+
  theme_bw(base_size = 14)+
  ylim(c(-3.5, 3.5))+
  ggtitle('week 4')+
  scale_color_manual(values=c("#5E3C99", 'purple', "#FDB863", "red", 'grey'))+
  scale_fill_manual(values=c("#5E3C99", 'purple', "#FDB863", "red", 'grey'))+
  xlim(c(30, 85))


w6 <-autoplot(aaW6, nrow=1, xlab='age (days)', ylab=expression(paste("hazard (", italic("H"), ")")))+
  geom_abline(intercept = 0, slope=0, linetype='dashed')+
  theme_bw(base_size = 14)+
  ylim(c(-1.7, 1.7))+
  ggtitle('week 6')+
  scale_color_manual(values=c("#5E3C99", 'purple', "#FDB863", "red", 'grey'))+
  scale_fill_manual(values=c("#5E3C99", 'purple', "#FDB863", "red", 'grey'))+
  xlim(c(44, 85))

ggarrange(w4, w6, ncol=1, common.legend = T, legend='right')


 
