library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(ggpubr)

masterFile='/Users/ar_ni/OneDrive/Desktop/APOE_MWM/apoe22_33_44_mwm_combined.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/APOE_MWM/R_Graphs/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df2<-subset(df, (APOE=='APOE2/2')) #Keeps only Genotypes 2/2, 3/3, and 4/4
df3<-subset(df, (APOE=='APOE3/3'))
df4<-subset(df, (APOE=='APOE4/4'))
dfAll<-rbind(df2, df3, df4)

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
dfNW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$NW.time, dfProbe$NW.distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW, quadrant='NW')
colnames(dfNW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$SW.time, dfProbe$SW.distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW, quadrant='SW')
colnames(dfSW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$NE.time, dfProbe$NE.distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE, quadrant='NE')
colnames(dfNE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$SE.time, dfProbe$SE.distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE, quadrant='SE')
colnames(dfSE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')

dfQuad<-rbind(dfSW, dfSE, dfNW, dfNE)
dfQuad2<-subset(dfQuad, (Genotype=='APOE2/2')) #Info on all quadrants for genotype APOE2/2
dfQuad3<-subset(dfQuad, (Genotype=='APOE3/3')) #Info on all quadrants for genotype APOE3/3
dfQuad4<-subset(dfQuad, (Genotype=='APOE4/4')) #Info on all quadrants for genotype APOE4/4
dfQuad2<-na.omit(dfQuad2) #Omit n/a
dfQuad3<-na.omit(dfQuad3) #Omit n/a
dfQuad4<-na.omit(dfQuad4) #Omit n/a
dfQuadp1<-subset(dfQuad, Day=='Probe_D5')
dfQuadp2<-subset(dfQuad, Day=='Probe_D8')

#Adjust SW Time for covariate Mean Speed
cor.test(dfAveraged$SW.time, dfAveraged$Mean.speed)
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
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
ggsave(paste(outpath,'poolTimeAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


data.lm <- lm(Duration ~ APOE + Stage, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Duration ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolTimeAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfAveraged, Stage %in% "Day4")

data.lm <- lm(Duration ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Duration ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolTimeAPOEDay4stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#__________________________________________________________________

#2. Time to platform over acquisition day for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='Duration', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top',
       facet.by="Sex"
)
ggsave(paste(outpath,'poolTimeAPOESex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempFemales <- subset(dfAveraged, Sex %in% "F")
tempMales <- subset(dfAveraged, Sex %in% "M")

dataF.lm <- lm(Duration ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(Duration ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(Duration ~ APOE, data = tempFemales)
testMethodM<-oneway.test(Duration ~ APOE, data = tempMales)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$statistic, testMethodF$p.value, testMethodF$parameter[1], nrow(tempFemales)) #Get values from summary
)

mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$statistic, testMethodM$p.value, testMethodM$parameter[1], nrow(tempMales)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTableF<-matrix(nrow=5, ncol=5)
postHocTableF[1,]=c('', '', '', '', '')
postHocTableF[2,]=c('Female TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolTimeAPOESexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfAveraged, Stage %in% "Day2")
temp <- subset(temp, Sex %in% "M")

data.lm <- lm(Duration ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Duration ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(paste("Male", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolTimeAPOEDay2Malestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfAveraged, Stage %in% "Day2")
temp <- subset(temp, APOE %in% "APOE3/3")

testMethod<-t.test(Duration ~ Sex, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeAPOE33Day2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Stage %in% "Day5")
temp <- subset(temp, APOE %in% "APOE4/4")

testMethod<-t.test(Duration ~ Sex, data = temp)


mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeAPOE44Day5stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________________________

#3. Distance to Platform over acquisition day for APOE 2/2, 3/3, 4/4
ggline(dfAveraged, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

data.lm <- lm(Distance ~ APOE, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#__________________________________________________________________

#4. Distance to Platform over acquisition day for APOE 2/2, 3/3, 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top',
       facet.by="Sex"
)
ggsave(paste(outpath,'poolDistAPOESex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

tempFemales <- subset(dfAveraged, Sex %in% "F")
tempMales <- subset(dfAveraged, Sex %in% "M")

dataF.lm <- lm(Distance ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(Distance ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(Distance ~ APOE, data = tempFemales)
testMethodM<-oneway.test(Distance ~ APOE, data = tempMales)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$statistic, testMethodF$p.value, testMethodF$parameter[1], nrow(tempFemales)) #Get values from summary
)

mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$statistic, testMethodM$p.value, testMethodM$parameter[1], nrow(tempMales)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTableF<-matrix(nrow=5, ncol=5)
postHocTableF[1,]=c('', '', '', '', '')
postHocTableF[2,]=c('Female TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOESexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#__________________________________________________________________

#5. NORMALIZED distance swam in SW quadrant for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(NormSWDist ~ APOE, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(NormSWDist ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfAveraged, Stage %in% "Day4")

data.lm <- lm(NormSWDist ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(NormSWDist ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEDay4stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#___________________________________________________________________

#6. NORMALIZED Time swam in SW quadrant for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='NormSWTime', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top')
ggsave(paste(outpath,'NormSWTimeAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(NormSWTime ~ APOE, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(NormSWTime ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'NormSWTimeAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#_________________________________________________________________________

#7. Mean speed for APOE 2/2, 3/3, and 4/4
ggline(dfAveraged, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top')
ggsave(paste(outpath,'MeanSpeedAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(Mean.speed ~ APOE, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Mean.speed ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'MeanSpeedAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <-subset(dfAveraged, Stage %in% "Day1")

data.lm <- lm(NormSWTime ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Mean.speed ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'MeanSpeedAPOEDay1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#_________________________________________________________________________
#8. SW Time swam adjusted for Mean Speed
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'SWTimeAdjustAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(residuals ~ APOE, data = dfAveraged)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(residuals ~ APOE, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'SWTimeAdjustAPOEstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfAveraged, Stage %in% "Day5")

data.lm <- lm(residuals ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(residuals ~ APOE, data = temp)
postHoc<-pairwise.t.test(temp$residuals, temp$APOE)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'SWTimeAdjustAPOEDay5stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS
#9. Barplot with Standard Error of DISTANCE in Quadrant for Genotype 2/2-- PROBE TRIAL
myMean<-aggregate(dfQuad2$Distance, by=list(Day=dfQuad2$Day, Quadrant=dfQuad2$Quadrant), mean)
mySD<-aggregate(dfQuad2$Distance, by=list(Day=dfQuad2$Day, Quadrant=dfQuad2$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)
#myMean$SE<-mySD$x/sqrt(nrow(subset(dfQuad2, (dfQuad2$Quadrant==mySD$Quadrant))))

head(subset(dfQuad2, (dfQuad2$Quadrant==mySD$Quadrant)), 25)
head(myMean, 25)
head(mySD, 25)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuad_2APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant APOE2/2",
                    ylab="Distance (m)",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Quadrant",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE * 2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of DISTANCE in Quadrant for Genotype 3/3-- PROBE TRIAL
myMean<-aggregate(dfQuad3$Distance, by=list(Day=dfQuad3$Day, Quadrant=dfQuad3$Quadrant), mean)
mySD<-aggregate(dfQuad3$Distance, by=list(Day=dfQuad3$Day, Quadrant=dfQuad3$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuad_3APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant APOE3/3",
                    ylab="Distance (m)",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Quadrant",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE * 2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of DISTANCE in Quadrant for Genotype 4/4-- PROBE TRIAL
myMean<-aggregate(dfQuad4$Distance, by=list(Day=dfQuad4$Day, Quadrant=dfQuad4$Quadrant), mean)
mySD<-aggregate(dfQuad4$Distance, by=list(Day=dfQuad4$Day, Quadrant=dfQuad4$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuad_4APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant APOE4/4",
                    ylab="Distance (m)",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Quadrant",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE * 2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Distance in Quadrant SW for All Genotypes-- PROBE TRIAL
myMean<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Genotype=dfSW$Genotype), mean)
mySD<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Genotype=dfSW$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6:7,1:3]
tabbedSE<-tabbedSE[6:7,1:3]

pdf(file='ProbeDistInQuadSW.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance in Quadrant SW",
                    ylab="Distance (m)",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Day",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file
