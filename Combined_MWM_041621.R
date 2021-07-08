library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(ggpubr)

masterFile='/Users/ar_ni/Scripts/apoe22_33_44__2HN_3HN_4HN_mwm_combined.csv'
outpath='/Users/ar_ni/Scripts/MWM_Combined_Graphs/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)
df<-subset(df, (Sex!='FHet'))

df2<-subset(df, (APOE=='E22'))
df3<-subset(df, (APOE=='E33'))
df4<-subset(df, (APOE=='E44'))
df2hn<-subset(df, (APOE=='E2HN'))
df3hn<-subset(df, (APOE=='E3HN'))
df4hn<-subset(df, (APOE=='E4HN'))
dfAll<- rbind(df2, df3, df4, df2hn, df3hn, df4hn)

#Correct for inaccurate distance calibration in AnyMaze
dfAll$Distance=dfAll$Distance/10
dfAll$NE.distance=dfAll$NE.distance/10
dfAll$NW.distance=dfAll$NW.distance/10
dfAll$SE.distance=dfAll$SE.distance/10
dfAll$SW.distance=dfAll$SW.distance/10

#Remove Probe Trials from Data Set
df1<-subset(dfAll, (Stage=='Day1'))
df2<-subset(dfAll, (Stage=='Day2'))
df3<-subset(dfAll, (Stage=='Day3'))
df4<-subset(dfAll, (Stage=='Day4'))
df5<-subset(dfAll, (Stage=='Day5'))
dfFin<-rbind(df1,df2,df3,df4,df5) #dfFin contains info from all regular trials

#Normalize time and distance in target region
dfFin$NormSWTime<-dfFin$SW.time/dfFin$Duration
dfFin$NormSWDist<-dfFin$SW.distance/dfFin$Distance
dfFin<-subset(dfFin, (NormSWTime <= 1))

#Averages 4 trials per day for each mouse
dfAveraged<-aggregate(.~Animal+APOE+Sex+Stage, dfFin, mean, na.action=na.pass)
#dfAveraged<-aggregate(dfFin, list(dfFin$Animal,dfFin$APOE,dfFin$Sex,dfFin$Stage), mean)
dfAPOE<-subset(dfAveraged, (APOE=='E22' | APOE=='E33' | APOE=='E44'))
dfHN<-subset(dfAveraged, (APOE=='E2HN' | APOE=='E3HN' | APOE=='E4HN'))

#Extract probe trials only
dfp1<-subset(dfAll, (Stage=='Probe_D5'))
dfp2<-subset(dfAll, (Stage=='Probe_D8'))
dfProbe<-rbind(dfp1, dfp2) #dfProbe contains info from all probe trials

#Normalize Probe Distances
dfProbe$DistTot<-dfProbe$NE.distance+dfProbe$NW.distance+dfProbe$SE.distance+dfProbe$SW.distance
dfProbe$NE.Dist.Norm<-dfProbe$NE.distance/dfProbe$DistTot
dfProbe$NW.Dist.Norm<-dfProbe$NW.distance/dfProbe$DistTot
dfProbe$SE.Dist.Norm<-dfProbe$SE.distance/dfProbe$DistTot
dfProbe$SW.Dist.Norm<-dfProbe$SW.distance/dfProbe$DistTot

#Normalize Probe Times
dfProbe$TimeTot<-dfProbe$NE.time+dfProbe$NW.time+dfProbe$SE.time+dfProbe$SW.time
dfProbe$NE.Time.Norm<-dfProbe$NE.time/dfProbe$TimeTot
dfProbe$NW.Time.Norm<-dfProbe$NW.time/dfProbe$TimeTot
dfProbe$SE.Time.Norm<-dfProbe$SE.time/dfProbe$TimeTot
dfProbe$SW.Time.Norm<-dfProbe$SW.time/dfProbe$TimeTot

#This section adds separate information by quadrants and adds quadrant labels
dfNW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Diet, dfProbe$Stage, dfProbe$Sex, dfProbe$NW.time, dfProbe$NW.distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW, quadrant='NW')
colnames(dfNW)<-c('Animal', 'Genotype', 'Diet', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Diet, dfProbe$Stage, dfProbe$Sex, dfProbe$SW.time, dfProbe$SW.distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW, quadrant='SW')
colnames(dfSW)<-c('Animal', 'Genotype', 'Diet', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Diet, dfProbe$Stage, dfProbe$Sex, dfProbe$NE.time, dfProbe$NE.distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE, quadrant='NE')
colnames(dfNE)<-c('Animal', 'Genotype', 'Diet', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Diet, dfProbe$Stage, dfProbe$Sex, dfProbe$SE.time, dfProbe$SE.distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE, quadrant='SE')
colnames(dfSE)<-c('Animal', 'Genotype', 'Diet', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')

dfSWp1<-subset(dfSW, (Day=='Probe_D5'))
dfSWp2<-subset(dfSW, (Day=='Probe_D8'))

dfQuad<-rbind(dfSW, dfSE, dfNW, dfNE)
dfQuad2<-subset(dfQuad, (Genotype=='E22')) #Info on all quadrants for genotype APOE2/2
dfQuad3<-subset(dfQuad, (Genotype=='E33')) #Info on all quadrants for genotype APOE3/3
dfQuad4<-subset(dfQuad, (Genotype=='E44')) #Info on all quadrants for genotype APOE4/4
dfQuad2hn<-subset(dfQuad, (Genotype=='E2HN')) #Info on all quadrants for genotype E2HN
dfQuad3hn<-subset(dfQuad, (Genotype=='E3HN')) #Info on all quadrants for genotype E3HN
dfQuad4hn<-subset(dfQuad, (Genotype=='E4HN')) #Info on all quadrants for genotype E4HN
dfQuad2<-na.omit(dfQuad2) #Omit n/a
dfQuad3<-na.omit(dfQuad3) #Omit n/a
dfQuad4<-na.omit(dfQuad4) #Omit n/a
dfQuad2hn<-na.omit(dfQuad2hn) #Omit n/a
dfQuad3hn<-na.omit(dfQuad3hn) #Omit n/a
dfQuad4hn<-na.omit(dfQuad4hn) #Omit n/a
dfQuadp1<-subset(dfQuad, Day=='Probe_D5')
dfQuadp2<-subset(dfQuad, Day=='Probe_D8')

#Adjust Norm SW Time for covariate Mean Speed
cor.test(dfAveraged$NormSWTime, dfAveraged$Mean.speed)
lm1<-lm(SW.time ~ Mean.speed, data=dfAveraged)
cor.test(lm1$residuals, dfAveraged$Mean.speed)
dfAveraged<-cbind(dfAveraged, residuals=lm1$residuals)

#__________________________________________________________________________
#Rewrite Data frame for Stats
write.csv(dfAveraged, file='mwmStatsAvgAPOE.csv')
write.csv(dfQuadp1, file='mwmProbeTrial1APOE.csv')
write.csv(dfQuadp2, file='mwmProbeTrial2APOE.csv')

#START PLOTTING
#__________________________________________________________________________

#1. Time to platform over acquisition day for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTime.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

data.lm <- lm(Duration ~ APOE + Stage, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Duration ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
        cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=11, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[4,1], tukey.test$APOE[4,2], tukey.test$APOE[4,3], tukey.test$APOE[4,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[11,1], tukey.test$APOE[11,2], tukey.test$APOE[11,3], tukey.test$APOE[11,4])
postHocTable[6,]=c('E2HN-E3HN', tukey.test$APOE[7,1], tukey.test$APOE[7,2], tukey.test$APOE[7,3], tukey.test$APOE[7,4])
postHocTable[7,]=c('E4HN-E2HN', tukey.test$APOE[9,1], tukey.test$APOE[9,2], tukey.test$APOE[9,3], tukey.test$APOE[9,4])
postHocTable[8,]=c('E4HN-E3HN', tukey.test$APOE[12,1], tukey.test$APOE[12,2], tukey.test$APOE[12,3], tukey.test$APOE[12,4])
postHocTable[9,]=c('APOE2/2-E2HN', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[10,]=c('APOE3/3-E3HN', tukey.test$APOE[10,1], tukey.test$APOE[10,2], tukey.test$APOE[10,3], tukey.test$APOE[10,4])
postHocTable[11,]=c('APOE4/4-E4HN', tukey.test$APOE[15,1], tukey.test$APOE[15,2], tukey.test$APOE[15,3], tukey.test$APOE[15,4])

myfile<-paste(outpath,'poolTimestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeHN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeAPOE22_E2HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeAPOE33_E3HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeAPOE44_E4HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#2. Time to platform over acquisition day for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolTimeAPOESex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#2. Time to platform over acquisition day for APOE 2/2, 3/3, and 4/4, separated by diet
ggline(dfAveraged, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E2HN vs E3HN vs E4HN
ggline(dfHN, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeHNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#APOE22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeAPOE22_E2HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


#APOE33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeAPOE33_E3HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempdietC<-subset(temp, (Diet=='1'))
tempdietHFD<-subset(temp, (Diet=='2'))

tempstats<-subset(temp, (Diet=='1' & Stage=='Day2'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Diet=='2' & Stage=='Day2'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#APOE44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolTimeAPOE44_E4HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempdietC<-subset(temp, (Diet=='1'))
tempdietHFD<-subset(temp, (Diet=='2'))

tempstats<-subset(temp, (Diet=='1' & Stage=='Day5'))
data.lm <- lm(Duration ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]
#__________________________________________________________________

#2. Distance to Platform over acquisition day for APOE 2/2, 3/3, 4/4
ggline(dfAveraged, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette =c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDist.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

data.lm <- lm(Distance ~ APOE + Stage, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
        cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=11, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[4,1], tukey.test$APOE[4,2], tukey.test$APOE[4,3], tukey.test$APOE[4,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[11,1], tukey.test$APOE[11,2], tukey.test$APOE[11,3], tukey.test$APOE[11,4])
postHocTable[6,]=c('E2HN-E3HN', tukey.test$APOE[7,1], tukey.test$APOE[7,2], tukey.test$APOE[7,3], tukey.test$APOE[7,4])
postHocTable[7,]=c('E4HN-E2HN', tukey.test$APOE[9,1], tukey.test$APOE[9,2], tukey.test$APOE[9,3], tukey.test$APOE[9,4])
postHocTable[8,]=c('E4HN-E3HN', tukey.test$APOE[12,1], tukey.test$APOE[12,2], tukey.test$APOE[12,3], tukey.test$APOE[12,4])
postHocTable[9,]=c('APOE2/2-E2HN', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[10,]=c('APOE3/3-E3HN', tukey.test$APOE[10,1], tukey.test$APOE[10,2], tukey.test$APOE[10,3], tukey.test$APOE[10,4])
postHocTable[11,]=c('APOE4/4-E4HN', tukey.test$APOE[15,1], tukey.test$APOE[15,2], tukey.test$APOE[15,3], tukey.test$APOE[15,4])

myfile<-paste(outpath,'poolDiststats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E2HN vs E3HN vs E4HN
ggline(dfHN, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistHN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE22_E2HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE33_E3HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day1'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE44_E4HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day1'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#2. Distance to Platform over acquisition day for APOE 2/2, 3/3, 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette =c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolDistAPOESex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#2. Distance to Platform over acquisition day for APOE 2/2, 3/3, 4/4, separated by diet
ggline(dfAveraged, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette =c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#E2HN vs E3HN vs E4HN
ggline(dfHN, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistHNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#APOE22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistAPOE22_E2HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


#APOE33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistAPOE33_E3HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempdietC<-subset(temp, (Diet=='1'))
tempdietHFD<-subset(temp, (Diet=='2'))

tempstats<-subset(temp, (Diet=='1' & Stage=='Day3'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Diet=='2' & Stage=='Day1'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Diet=='2' & Stage=='Day2'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#APOE44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'poolDistAPOE44_E4HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempdietC<-subset(temp, (Diet=='1'))
tempdietHFD<-subset(temp, (Diet=='2'))

tempstats<-subset(temp, (Diet=='1' & Stage=='Day5'))
data.lm <- lm(Distance ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#__________________________________________________________________

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDist.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistHN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE22_E2HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE33_E3HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE44_E4HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWDistSex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, 3/3, and 4/4, separated by diet
ggline(dfAveraged, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistDiet.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistHNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistAPOE22_E2HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistAPOE33_E3HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWDistAPOE44_E4HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#___________________________________________________________________
#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='NormSWTime', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top')
ggsave(paste(outpath,'NormSWTimeAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='NormSWTime', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWTimeAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, 3/3, and 4/4, separated by diet
ggline(dfAveraged, x='Stage', y='NormSWTime', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAPOEDiet.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#_________________________________________________________________________

#7. Mean speed 
ggline(dfAveraged, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top')
ggsave(paste(outpath,'MeanSpeed.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'MeanSpeedAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'MeanSpeedHN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'MeanSpeedAPOE22_E2HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'MeanSpeedAPOE33_E3HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'MeanSpeedAPOE44_E4HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#7. Mean speed, separated by sex
ggline(dfAveraged, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'MeanSpeedAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#7. Mean speed, separated by diet
ggline(dfAveraged, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedAPOEDiet.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedHNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedAPOE22_E2HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedAPOE33_E3HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by='Diet')
ggsave(paste(outpath,'MeanSpeedAPOE44_E4HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#_________________________________________________________________________

#5. Norm SW Time swam adjusted for Mean Speed
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjust.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustHN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE22_E2HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE33_E3HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE44_E4HN.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


#5. Norm SW Time swam adjusted for Mean Speed, separated by sex
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWTimeAdjustSex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#5. Norm SW Time swam adjusted for Mean Speed, separated by diet
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustDiet.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#E22 vs E33 vs E44
ggline(dfAPOE, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustAPOEDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfAPOE, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfAPOE, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E2HN vs E3HNvs E4HN
ggline(dfHN, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustHNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(dfHN, (Stage=='Day1'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(dfHN, (Stage=='Day3'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E22 vs E2HN
temp<-subset(dfAveraged, (APOE=='E22' | APOE=='E2HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE22_E2HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E33 vs E3HN
temp<-subset(dfAveraged, (APOE=='E33' | APOE=='E3HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE33_E3HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempstats<-subset(temp, (Stage=='Day2'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

tempstats<-subset(temp, (Stage=='Day5'))
data.lm <- lm(Mean.speed ~ APOE, data = tempstats)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)
tukey.test$APOE[1,4]

#E44 vs E4HN
temp<-subset(dfAveraged, (APOE=='E44' | APOE=='E4HN'))
ggline(temp, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by='Diet')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE44_E4HNDiet.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#_________________________________________________________________________
#8. Winding Numbers for each genotype
ggline(dfAveraged, x='Stage', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top')
ggsave(paste(outpath,'WindingAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#8. Winding Numbers for each genotype, separated by sex
ggline(dfAveraged, x='Stage', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top', facet.by='Sex')
ggsave(paste(outpath,'WindingAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#8. Winding Numbers for each genotype, separated by diet
ggline(dfAveraged, x='Stage', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top', facet.by='Diet')
ggsave(paste(outpath,'WindingAPOEDiet.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)
#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS
#Distance in each Quadrant for APOE22 in Probes
ggline(dfQuad2, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE22Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in each Quadrant for APOE33 in Probes
ggline(dfQuad3, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE33Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in each Quadrant for APOE44 in Probes
ggline(dfQuad4, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE44Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in each Quadrant for APOE2HN in Probes
ggline(dfQuad2hn, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE2HNProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in each Quadrant for APOE3HN in Probes
ggline(dfQuad3hn, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE3HNProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in each Quadrant for APOE4HN in Probes
ggline(dfQuad4hn, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE4HNProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Distance in SW for all genotypes in Probes
ggline(dfSW, x='Genotype', y='Distance', color='Genotype', fill='Genotype',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'SWDistAPOEProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Normalized Distance in SW for all genotypes in Probes
ggline(dfSW, x='Genotype', y='DistNorm', color='Genotype', fill='Genotype',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red','pink','gray','blue'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'NormSWDistAPOEProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

