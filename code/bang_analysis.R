library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggpubr)
library(smatr)


rm(list=ls())
getwd()
setwd("/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/Claudia Sun/Tom Nonacs/bang assay_TomNonacs")

dir()
d <- read.csv("3-31_climbing_seperation.csv")

head(d)
range(d$n)

table(d$climbing_group)


d$Climbing <-  as.factor(d$climbing_group)
levels(d$Climbing) <- c('bottom', 'mid-bottom', 'mid-bottom', 'mid-top', 'mid-bottom', 'mid-top', 'mid-top', 'top') 
d$ClimbingFraction <- d$Climbing
levels(d$Climbing)
levels(d$ClimbingFraction) <- c('0/3', '1/3', '2/3', '3/3')

d$ClimbingFraction <- factor(d$ClimbingFraction, levels=rev(levels(d$ClimbingFraction))) 

table(d$Climbing, d$treatment)

d$prop_surv <- 1-d$prop_dead

ggplot(d, aes(y=prop_surv, x=Climbing, color=treatment, fill=treatment))+
  geom_boxplot(alpha=0.4, outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1))+
  ylab('proportion surviving')+
  theme_bw(base_size = 14)+
  ggtitle('bang sensitivity')

ggplot(d, aes(y=prop_surv, x=ClimbingFraction, color=treatment, fill=treatment))+
  geom_boxplot(alpha=0.4, outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1))+
  ylab('proportion surviving')+
  theme_bw(base_size = 14)+
  ggtitle('bang sensitivity')

d$ClimbingFraction <- factor(d$ClimbingFraction, levels=rev(levels(d$ClimbingFraction)))
d$ClimbingFraction

table(d$ClimbingFraction, d$Climbing)

bangPlot <- ggplot(d, aes(y=prop_surv, x=ClimbingFraction, color=Climbing, fill=Climbing))+
  geom_boxplot(alpha=0.4, outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1))+
  ylab('proportion surviving')+
  theme_bw(base_size = 14)+
  facet_wrap(~d$treatment )+
  xlab('')+
  scale_color_manual(values=rev(c("red", "#FDB863", 'purple', "#5E3C99"))) +
  scale_fill_manual(values=rev(c("red", "#FDB863", 'purple', "#5E3C99"))) 

bangPlot

save(d, bangPlot, file='bangplot')


aggregate(prop_surv ~ Climbing * treatment, d, mean)

table(d$climbing_group, d$treatment)

library(betareg)
library(lmtest)
?betareg

d$prop_surv

y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs} # function to transform 0-inflated proportional data (Smithson and Verkuilen 2006).  see: https://stackoverflow.com/questions/26385617/proportion-modeling-betareg-errors

par(mfrow=c(2,1))
hist(d$prop_dead)
hist(y.transf.betareg(d$prop_dead))

model.beta = betareg(y.transf.betareg(d$prop_dead) ~ Climbing * treatment, d)

m1 <- betareg(y.transf.betareg(d$prop_dead) ~ 1, d)
m2 <- betareg(y.transf.betareg(d$prop_dead) ~ Climbing, d)
m3 <- betareg(y.transf.betareg(d$prop_dead) ~ treatment, d)
m4 <- betareg(y.transf.betareg(d$prop_dead) ~ Climbing + treatment, d)
m5 <- betareg(y.transf.betareg(d$prop_dead) ~ Climbing * treatment, d)

lrtest(model.beta)

lrtest(m1, m2, m3, m4, m5) # the last line shows that model 5, the addition of an interaction term, is significant
lrtest(m1, m5)
lrtest(m4, m5)


