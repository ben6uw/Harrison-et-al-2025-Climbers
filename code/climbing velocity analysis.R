library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggpubr)
library(smatr)

rm(list=ls())

setwd('~/Google Drive/My Drive/Documents/Claudia Sun/Tom Nonacs/smurfs and climbing velocity')
dir()

v <- readxl::read_excel('mm_SmurfClimbing_Velocity_All_Days.xlsx')

head(v)
colnames(v)

"(velocity)_kalman_y_slope_(mm/frame)"

v$velocity <- v$`(velocity)_kalman_y_slope_(mm/frame)` * 30 # per frame distance in mm, frame rate was 30/s, so multiply these by 30 to get velocity in mm/s
hist(v$velocity)

#[NOTE]: Flies were originally cleared on 12/27/2023 at around 8pm, Collected on 12/30/2023 at around 8pm, and sexed on 1/2/2024.  So, 8am on 12/29 will serve as age = 0days
start <- "2023-12-29 UTC"
v$age <- as.numeric(difftime(as.POSIXct( v$day, format="%Y-%m-%dT%H:%M"), as.POSIXct(start), units="day"))
table(v$age)

head(v)
v$group <- as.factor(v$group)
levels(v$group) <- c('bottom', 'top')

ggplot(v, aes(y=velocity, x=age, color=group, fill=group, group=age))+
  geom_jitter(width=0.1, size=0.1)+
  theme_bw(base_size = 16)+
  stat_summary(aes(group=group), fun.y=mean, geom="line", colour="#5E3C99")

head(v)
means <- v %>% group_by(group, age) %>% summarise(mean.velocity=mean(velocity), sd.velocity=sd(velocity), n = n()) 
head(means)
means$se.velocity <- means$sd.velocity/sqrt(means$n)
head(means)

means$fraction <- means$group
means$fraction <- factor(means$fraction, levels=c('top', 'bottom'))

ggplot(means, aes(y=mean.velocity, x=age, color=fraction , group=fraction))+
  geom_errorbar(aes(ymin=mean.velocity - se.velocity, ymax=mean.velocity + se.velocity), width=0.2, col='grey40')+
  geom_point()+
  geom_line()+
  scale_color_manual(values=rev(c("#5E3C99", 'red')))+
  ylab('mean velocity (mm/s) +/-se')+
  xlab('age (days)')+
  ylim(0, 11)+
  labs(color='climbing')+
  theme_bw(base_size = 16)


range(means$n)
means
str(v)

head(v)
table(v$pos)
table(v$batch)
table(v$vial)

summary(lmer(velocity ~ group * as.factor(age) + pos + batch + (1|vial), v))




