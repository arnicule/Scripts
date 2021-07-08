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
dfNW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Winding, dfProbe$NW.time, dfProbe$NW.distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW, quadrant='NW')
colnames(dfNW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Winding', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Winding, dfProbe$SW.time, dfProbe$SW.distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW, quadrant='SW')
colnames(dfSW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Winding', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Winding, dfProbe$NE.time, dfProbe$NE.distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE, quadrant='NE')
colnames(dfNE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Winding', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Winding, dfProbe$SE.time, dfProbe$SE.distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE, quadrant='SE')
colnames(dfSE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Winding', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')

dfSWp1<-subset(dfSW, (Day=='Probe_D5'))
dfSWp2<-subset(dfSW, (Day=='Probe_D8'))

dfQuad<-rbind(dfSW, dfSE, dfNW, dfNE)
dfQuad2<-subset(dfQuad, (Genotype=='APOE2/2')) #Info on all quadrants for genotype APOE2/2
dfQuad3<-subset(dfQuad, (Genotype=='APOE3/3')) #Info on all quadrants for genotype APOE3/3
dfQuad4<-subset(dfQuad, (Genotype=='APOE4/4')) #Info on all quadrants for genotype APOE4/4
dfQuad2<-na.omit(dfQuad2) #Omit n/a
dfQuad3<-na.omit(dfQuad3) #Omit n/a
dfQuad4<-na.omit(dfQuad4) #Omit n/a
dfQuadp1<-subset(dfQuad, Day=='Probe_D5')
dfQuadp2<-subset(dfQuad, Day=='Probe_D8')

#Adjust Normalized SW Time for covariate Mean Speed
cor.test(dfAveraged$NormSWTime, dfAveraged$Mean.speed)
lm1<-lm(NormSWTime ~ Mean.speed, data=dfAveraged)
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

#5. NORMALIZED distance swam in SW quadrant for APOE 2/2, 3/3, and 4/4 separated by sex
ggline(dfAveraged, x='Stage', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWDistAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

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

#6. NORMALIZED Time swam in SW quadrant for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='NormSWTime', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWTimeAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

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

#7. Mean speed for APOE 2/2, 3/3, and 4/4, separated by sex
ggline(dfAveraged, x='Stage', y='Mean.speed', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'MeanSpeedAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfAveraged, Sex %in% "F")
tempMales <- subset(dfAveraged, Sex %in% "M")

dataF.lm <- lm(Mean.speed ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(Mean.speed ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(Mean.speed ~ APOE, data = tempFemales)
testMethodM<-oneway.test(Mean.speed ~ APOE, data = tempMales)

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'meanSpeedAPOESexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)


temp <- subset(dfAveraged, Stage %in% "Day1")
temp <- subset(temp, Sex %in% "F")

data.lm <- lm(Mean.speed ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Mean.speed ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Female TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'meanSpeedAPOEDay1Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
#_________________________________________________________________________
#8. Normalized SW Time swam adjusted for Mean Speed
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
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

#8. Normalized SW Time swam adjusted for Mean Speed, separated by sex
ggline(dfAveraged, x='Stage', y='residuals', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWTimeAdjustAPOESex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#8. Winding Numbers for each genotype
ggline(dfAveraged, x='Stage', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top')
ggsave(paste(outpath,'WindingAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#8. Winding Numbers for each genotype, separated by sex
ggline(dfAveraged, x='Stage', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top', facet.by='Sex')
ggsave(paste(outpath,'WindingAPOEsex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

temp <- subset(dfAveraged, Stage %in% "Day5")
temp <- subset(temp, Sex %in% "F")

data.lm <- lm(Winding ~ APOE, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Winding ~ APOE, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Female TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'windingAPOEDay5Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS

#Distance in each Quadrant for APOE22 in Probes
ggline(dfQuad2, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE22Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#9. Barplot with Standard Error of DISTANCE in Quadrant for Genotype 2/2-- PROBE TRIAL
myMean<-aggregate(dfQuad2$Distance, by=list(Day=dfQuad2$Day, Quadrant=dfQuad2$Quadrant), mean)
mySD<-aggregate(dfQuad2$Distance, by=list(Day=dfQuad2$Day, Quadrant=dfQuad2$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

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

dfQuad2p1<-subset(dfQuadp1, Genotype=='APOE2/2')
dfQuad2p2<-subset(dfQuadp2, Genotype=='APOE2/2')

datap1.lm <- lm(Distance ~ Quadrant, data = dfQuad2p1)
datap1.aov <- aov(datap1.lm)
tukeyp1.test <- TukeyHSD(datap1.aov)
tukeyp1.test

datap2.lm <- lm(Distance ~ Quadrant, data = dfQuad2p2)
datap2.aov <- aov(datap2.lm)
tukeyp2.test <- TukeyHSD(datap2.aov)

testMethodp1<-oneway.test(Distance ~ Quadrant, data = dfQuad2p1)
testMethodp2<-oneway.test(Distance ~ Quadrant, data = dfQuad2p2)

mytTablep1<-as_tibble(
  cbind(paste("Probe 1", testMethodp1$data.name, sep=" "), testMethodp1$statistic, testMethodp1$p.value, testMethodp1$parameter[1], nrow(dfQuad2p1)) #Get values from summary
)

mytTablep2<-as_tibble(
  cbind(paste("Probe 2", testMethodp2$data.name, sep=" "), testMethodp2$statistic, testMethodp2$p.value, testMethodp2$parameter[1], nrow(dfQuad2p2)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTablep1<-matrix(nrow=8, ncol=5)
postHocTablep1[1,]=c('', '', '', '', '')
postHocTablep1[2,]=c('Probe 1 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep1[3,]=c('SE-SW', tukeyp1.test$Quadrant[1,1], tukeyp1.test$Quadrant[1,2], tukeyp1.test$Quadrant[1,3], tukeyp1.test$Quadrant[1,4])
postHocTablep1[4,]=c('NW-SW', tukeyp1.test$Quadrant[2,1], tukeyp1.test$Quadrant[2,2], tukeyp1.test$Quadrant[2,3], tukeyp1.test$Quadrant[2,4])
postHocTablep1[5,]=c('NE-SW', tukeyp1.test$Quadrant[3,1], tukeyp1.test$Quadrant[3,2], tukeyp1.test$Quadrant[3,3], tukeyp1.test$Quadrant[3,4])
postHocTablep1[6,]=c('NW-SE', tukeyp1.test$Quadrant[4,1], tukeyp1.test$Quadrant[4,2], tukeyp1.test$Quadrant[4,3], tukeyp1.test$Quadrant[4,4])
postHocTablep1[7,]=c('NE-SE', tukeyp1.test$Quadrant[5,1], tukeyp1.test$Quadrant[5,2], tukeyp1.test$Quadrant[5,3], tukeyp1.test$Quadrant[5,4])
postHocTablep1[8,]=c('NE-NW', tukeyp1.test$Quadrant[6,1], tukeyp1.test$Quadrant[6,2], tukeyp1.test$Quadrant[6,3], tukeyp1.test$Quadrant[6,4])

postHocTablep2<-matrix(nrow=8, ncol=5)
postHocTablep2[1,]=c('', '', '', '', '')
postHocTablep2[2,]=c('Probe 2 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep2[3,]=c('SE-SW', tukeyp2.test$Quadrant[1,1], tukeyp2.test$Quadrant[1,2], tukeyp2.test$Quadrant[1,3], tukeyp2.test$Quadrant[1,4])
postHocTablep2[4,]=c('NW-SW', tukeyp2.test$Quadrant[2,1], tukeyp2.test$Quadrant[2,2], tukeyp2.test$Quadrant[2,3], tukeyp2.test$Quadrant[2,4])
postHocTablep2[5,]=c('NE-SW', tukeyp2.test$Quadrant[3,1], tukeyp2.test$Quadrant[3,2], tukeyp2.test$Quadrant[3,3], tukeyp2.test$Quadrant[3,4])
postHocTablep2[6,]=c('NW-SE', tukeyp2.test$Quadrant[4,1], tukeyp2.test$Quadrant[4,2], tukeyp2.test$Quadrant[4,3], tukeyp2.test$Quadrant[4,4])
postHocTablep2[7,]=c('NE-SE', tukeyp2.test$Quadrant[5,1], tukeyp2.test$Quadrant[5,2], tukeyp2.test$Quadrant[5,3], tukeyp2.test$Quadrant[5,4])
postHocTablep2[8,]=c('NE-NW', tukeyp2.test$Quadrant[6,1], tukeyp2.test$Quadrant[6,2], tukeyp2.test$Quadrant[6,3], tukeyp2.test$Quadrant[6,4])

myfile<-paste(outpath,'ProbeQuadDistAPOE22stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Distance in each Quadrant for APOE33 in Probes
ggline(dfQuad3, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE33Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

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

dfQuad3p1<-subset(dfQuadp1, Genotype=='APOE3/3')
dfQuad3p2<-subset(dfQuadp2, Genotype=='APOE3/3')

datap1.lm <- lm(Distance ~ Quadrant, data = dfQuad3p1)
datap1.aov <- aov(datap1.lm)
tukeyp1.test <- TukeyHSD(datap1.aov)
tukeyp1.test

datap2.lm <- lm(Distance ~ Quadrant, data = dfQuad3p2)
datap2.aov <- aov(datap2.lm)
tukeyp2.test <- TukeyHSD(datap2.aov)

testMethodp1<-oneway.test(Distance ~ Quadrant, data = dfQuad3p1)
testMethodp2<-oneway.test(Distance ~ Quadrant, data = dfQuad3p2)

mytTablep1<-as_tibble(
  cbind(paste("Probe 1", testMethodp1$data.name, sep=" "), testMethodp1$statistic, testMethodp1$p.value, testMethodp1$parameter[1], nrow(dfQuad2p1)) #Get values from summary
)

mytTablep2<-as_tibble(
  cbind(paste("Probe 2", testMethodp2$data.name, sep=" "), testMethodp2$statistic, testMethodp2$p.value, testMethodp2$parameter[1], nrow(dfQuad2p2)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTablep1<-matrix(nrow=8, ncol=5)
postHocTablep1[1,]=c('', '', '', '', '')
postHocTablep1[2,]=c('Probe 1 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep1[3,]=c('SE-SW', tukeyp1.test$Quadrant[1,1], tukeyp1.test$Quadrant[1,2], tukeyp1.test$Quadrant[1,3], tukeyp1.test$Quadrant[1,4])
postHocTablep1[4,]=c('NW-SW', tukeyp1.test$Quadrant[2,1], tukeyp1.test$Quadrant[2,2], tukeyp1.test$Quadrant[2,3], tukeyp1.test$Quadrant[2,4])
postHocTablep1[5,]=c('NE-SW', tukeyp1.test$Quadrant[3,1], tukeyp1.test$Quadrant[3,2], tukeyp1.test$Quadrant[3,3], tukeyp1.test$Quadrant[3,4])
postHocTablep1[6,]=c('NW-SE', tukeyp1.test$Quadrant[4,1], tukeyp1.test$Quadrant[4,2], tukeyp1.test$Quadrant[4,3], tukeyp1.test$Quadrant[4,4])
postHocTablep1[7,]=c('NE-SE', tukeyp1.test$Quadrant[5,1], tukeyp1.test$Quadrant[5,2], tukeyp1.test$Quadrant[5,3], tukeyp1.test$Quadrant[5,4])
postHocTablep1[8,]=c('NE-NW', tukeyp1.test$Quadrant[6,1], tukeyp1.test$Quadrant[6,2], tukeyp1.test$Quadrant[6,3], tukeyp1.test$Quadrant[6,4])

postHocTablep2<-matrix(nrow=8, ncol=5)
postHocTablep2[1,]=c('', '', '', '', '')
postHocTablep2[2,]=c('Probe 2 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep2[3,]=c('SE-SW', tukeyp2.test$Quadrant[1,1], tukeyp2.test$Quadrant[1,2], tukeyp2.test$Quadrant[1,3], tukeyp2.test$Quadrant[1,4])
postHocTablep2[4,]=c('NW-SW', tukeyp2.test$Quadrant[2,1], tukeyp2.test$Quadrant[2,2], tukeyp2.test$Quadrant[2,3], tukeyp2.test$Quadrant[2,4])
postHocTablep2[5,]=c('NE-SW', tukeyp2.test$Quadrant[3,1], tukeyp2.test$Quadrant[3,2], tukeyp2.test$Quadrant[3,3], tukeyp2.test$Quadrant[3,4])
postHocTablep2[6,]=c('NW-SE', tukeyp2.test$Quadrant[4,1], tukeyp2.test$Quadrant[4,2], tukeyp2.test$Quadrant[4,3], tukeyp2.test$Quadrant[4,4])
postHocTablep2[7,]=c('NE-SE', tukeyp2.test$Quadrant[5,1], tukeyp2.test$Quadrant[5,2], tukeyp2.test$Quadrant[5,3], tukeyp2.test$Quadrant[5,4])
postHocTablep2[8,]=c('NE-NW', tukeyp2.test$Quadrant[6,1], tukeyp2.test$Quadrant[6,2], tukeyp2.test$Quadrant[6,3], tukeyp2.test$Quadrant[6,4])

myfile<-paste(outpath,'ProbeQuadDistAPOE33stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Distance in each Quadrant for APOE44 in Probes
ggline(dfQuad4, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'quadDistAPOE44Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

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

dfQuad4p1<-subset(dfQuadp1, Genotype=='APOE4/4')
dfQuad4p2<-subset(dfQuadp2, Genotype=='APOE4/4')

datap1.lm <- lm(Distance ~ Quadrant, data = dfQuad4p1)
datap1.aov <- aov(datap1.lm)
tukeyp1.test <- TukeyHSD(datap1.aov)
tukeyp1.test

datap2.lm <- lm(Distance ~ Quadrant, data = dfQuad4p2)
datap2.aov <- aov(datap2.lm)
tukeyp2.test <- TukeyHSD(datap2.aov)

testMethodp1<-oneway.test(Distance ~ Quadrant, data = dfQuad4p1)
testMethodp2<-oneway.test(Distance ~ Quadrant, data = dfQuad4p2)

mytTablep1<-as_tibble(
  cbind(paste("Probe 1", testMethodp1$data.name, sep=" "), testMethodp1$statistic, testMethodp1$p.value, testMethodp1$parameter[1], nrow(dfQuad2p1)) #Get values from summary
)

mytTablep2<-as_tibble(
  cbind(paste("Probe 2", testMethodp2$data.name, sep=" "), testMethodp2$statistic, testMethodp2$p.value, testMethodp2$parameter[1], nrow(dfQuad2p2)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTablep1<-matrix(nrow=8, ncol=5)
postHocTablep1[1,]=c('', '', '', '', '')
postHocTablep1[2,]=c('Probe 1 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep1[3,]=c('SE-SW', tukeyp1.test$Quadrant[1,1], tukeyp1.test$Quadrant[1,2], tukeyp1.test$Quadrant[1,3], tukeyp1.test$Quadrant[1,4])
postHocTablep1[4,]=c('NW-SW', tukeyp1.test$Quadrant[2,1], tukeyp1.test$Quadrant[2,2], tukeyp1.test$Quadrant[2,3], tukeyp1.test$Quadrant[2,4])
postHocTablep1[5,]=c('NE-SW', tukeyp1.test$Quadrant[3,1], tukeyp1.test$Quadrant[3,2], tukeyp1.test$Quadrant[3,3], tukeyp1.test$Quadrant[3,4])
postHocTablep1[6,]=c('NW-SE', tukeyp1.test$Quadrant[4,1], tukeyp1.test$Quadrant[4,2], tukeyp1.test$Quadrant[4,3], tukeyp1.test$Quadrant[4,4])
postHocTablep1[7,]=c('NE-SE', tukeyp1.test$Quadrant[5,1], tukeyp1.test$Quadrant[5,2], tukeyp1.test$Quadrant[5,3], tukeyp1.test$Quadrant[5,4])
postHocTablep1[8,]=c('NE-NW', tukeyp1.test$Quadrant[6,1], tukeyp1.test$Quadrant[6,2], tukeyp1.test$Quadrant[6,3], tukeyp1.test$Quadrant[6,4])

postHocTablep2<-matrix(nrow=8, ncol=5)
postHocTablep2[1,]=c('', '', '', '', '')
postHocTablep2[2,]=c('Probe 2 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTablep2[3,]=c('SE-SW', tukeyp2.test$Quadrant[1,1], tukeyp2.test$Quadrant[1,2], tukeyp2.test$Quadrant[1,3], tukeyp2.test$Quadrant[1,4])
postHocTablep2[4,]=c('NW-SW', tukeyp2.test$Quadrant[2,1], tukeyp2.test$Quadrant[2,2], tukeyp2.test$Quadrant[2,3], tukeyp2.test$Quadrant[2,4])
postHocTablep2[5,]=c('NE-SW', tukeyp2.test$Quadrant[3,1], tukeyp2.test$Quadrant[3,2], tukeyp2.test$Quadrant[3,3], tukeyp2.test$Quadrant[3,4])
postHocTablep2[6,]=c('NW-SE', tukeyp2.test$Quadrant[4,1], tukeyp2.test$Quadrant[4,2], tukeyp2.test$Quadrant[4,3], tukeyp2.test$Quadrant[4,4])
postHocTablep2[7,]=c('NE-SE', tukeyp2.test$Quadrant[5,1], tukeyp2.test$Quadrant[5,2], tukeyp2.test$Quadrant[5,3], tukeyp2.test$Quadrant[5,4])
postHocTablep2[8,]=c('NE-NW', tukeyp2.test$Quadrant[6,1], tukeyp2.test$Quadrant[6,2], tukeyp2.test$Quadrant[6,3], tukeyp2.test$Quadrant[6,4])

myfile<-paste(outpath,'ProbeQuadDistAPOE44stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Distance in SW for all genotypes in Probes
ggline(dfSW, x='Genotype', y='Distance', color='Genotype', fill='Genotype',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Distance (m)', legend='top', facet.by='Day')
ggsave(paste(outpath,'SWDistAPOEProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Barplot with Standard Error of Distance in Quadrant SW for All Genotypes-- PROBE TRIAL
myMean<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Genotype=dfSW$Genotype), mean)
mySD<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Genotype=dfSW$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(58)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6:7,1:3]
tabbedSE<-tabbedSE[6:7,1:3]

pdf(file='ProbeSWDistAPOE.pdf')
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

dfQuad2SW<-subset(dfSW, Genotype=='APOE2/2')
dfQuad3SW<-subset(dfSW, Genotype=='APOE3/3')
dfQuad4SW<-subset(dfSW, Genotype=='APOE4/4')

data2.lm <- lm(Distance ~ Day, data = dfQuad2SW)
data2.aov <- aov(data2.lm)
tukey2.test <- TukeyHSD(data2.aov)

data3.lm <- lm(Distance ~ Day, data = dfQuad3SW)
data3.aov <- aov(data3.lm)
tukey3.test <- TukeyHSD(data3.aov)

data4.lm <- lm(Distance ~ Day, data = dfQuad4SW)
data4.aov <- aov(data4.lm)
tukey4.test <- TukeyHSD(data4.aov)

testMethod2<-oneway.test(Distance ~ Day, data = dfQuad2SW)
testMethod3<-oneway.test(Distance ~ Day, data = dfQuad3SW)
testMethod4<-oneway.test(Distance ~ Day, data = dfQuad4SW)

mytTable2<-as_tibble(
  cbind(paste("APOE 2/2", testMethod2$data.name, sep=" "), testMethod2$statistic, testMethod2$p.value, testMethod2$parameter[1], nrow(dfQuad2SW)) #Get values from summary
)
mytTable3<-as_tibble(
  cbind(paste("APOE 3/3", testMethod3$data.name, sep=" "), testMethod3$statistic, testMethod3$p.value, testMethod3$parameter[1], nrow(dfQuad3SW)) #Get values from summary
)
mytTable4<-as_tibble(
  cbind(paste("APOE 4/4", testMethod4$data.name, sep=" "), testMethod4$statistic, testMethod4$p.value, testMethod4$parameter[1], nrow(dfQuad4SW)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable2<-matrix(nrow=3, ncol=5)
postHocTable2[1,]=c('', '', '', '', '')
postHocTable2[2,]=c('APOE 2/2 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable2[3,]=c('Probe2-Probe1', tukey2.test$Day[1,1], tukey2.test$Day[1,2], tukey2.test$Day[1,3], tukey2.test$Day[1,4])

postHocTable3<-matrix(nrow=3, ncol=5)
postHocTable3[1,]=c('', '', '', '', '')
postHocTable3[2,]=c('APOE 3/3 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable3[3,]=c('Probe2-Probe1', tukey3.test$Day[1,1], tukey3.test$Day[1,2], tukey3.test$Day[1,3], tukey3.test$Day[1,4])

postHocTable4<-matrix(nrow=3, ncol=5)
postHocTable4[1,]=c('', '', '', '', '')
postHocTable4[2,]=c('APOE 4/4 TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable4[3,]=c('Probe2-Probe1', tukey4.test$Day[1,1], tukey4.test$Day[1,2], tukey4.test$Day[1,3], tukey4.test$Day[1,4])

myfile<-paste(outpath,'ProbeQuadDistAPOEstats.csv')
write.table(mytTable2, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTable3, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(mytTable4, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTable2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTable3, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTable4, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of Normalized Distance in Quadrant SW in Probe Day 1 for All Genotypes-- PROBE TRIAL
myMean<-aggregate(dfSWp1$DistNorm, by=list(Genotype=dfSWp1$Genotype), mean)
mySD<-aggregate(dfSWp1$DistNorm, by=list(Genotype=dfSWp1$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(29)

tabbedMeans<-tapply(myMean$x, list(myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6,1:3]
tabbedSE<-tabbedSE[6,1:3]

pdf(file='ProbeNormSWDistP1APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 1",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Genotype",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW in Probe Day 2 for All Genotypes-- PROBE TRIAL
myMean<-aggregate(dfSWp2$DistNorm, by=list(Genotype=dfSWp2$Genotype), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Genotype=dfSWp2$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(29)

tabbedMeans<-tapply(myMean$x, list(myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[7,1:3]
tabbedSE<-tabbedSE[7,1:3]

pdf(file='ProbeNormSWDistP2APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 2",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Genotype",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW in Probe Day 1 for All Genotypes-- PROBE TRIAL, separated by sex
myMean<-aggregate(dfSWp1$DistNorm, by=list(Sex=dfSWp1$Sex, Genotype=dfSWp1$Genotype), mean)
mySD<-aggregate(dfSWp1$DistNorm, by=list(Sex=dfSWp1$Sex, Genotype=dfSWp1$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(29)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[1:2,1:3]
tabbedSE<-tabbedSE[1:2,1:3]

pdf(file='ProbeNormSWDistP1APOESex.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 1",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Sex",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW in Probe Day 2 for All Genotypes-- PROBE TRIAL, separated by sex
myMean<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Genotype=dfSWp2$Genotype), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Genotype=dfSWp2$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(29)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[1:2,1:3]
tabbedSE<-tabbedSE[1:2,1:3]

pdf(file='ProbeNormSWDistP2APOESex.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 2",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Sex",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Winding Numbers in Probe Day 1 for APOE22/33/44
ggline(dfp1, x='APOE', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top')
ggsave(paste(outpath,'WindingAPOEProbe1.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(Winding ~ APOE, data = dfp1)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Winding ~ APOE, data = dfp1)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp1))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'WindingAPOEProbe1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Winding Numbers in Probe Day 1 for APOE22/33/44, separated by sex
ggline(dfp1, x='APOE', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top', facet.by='Sex')
ggsave(paste(outpath,'WindingAPOEProbe1sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp1, Sex %in% "F")
tempMales <- subset(dfp1, Sex %in% "M")

dataF.lm <- lm(Winding ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(Winding ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(Winding ~ APOE, data = tempFemales)
testMethodM<-oneway.test(Winding ~ APOE, data = tempMales)

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'windingAPOEProbe1sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Winding Numbers in Probe Day 2 for APOE22/33/44
ggline(dfp2, x='APOE', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top')
ggsave(paste(outpath,'WindingAPOEProbe2.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(Winding ~ APOE, data = dfp2)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Winding ~ APOE, data = dfp2)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp2))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'WindingAPOEProbe2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Winding Numbers in Probe Day 2 for APOE22/33/44, separated by sex
ggline(dfp2, x='APOE', y='Winding', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Winding Numbers', legend='top', facet.by='Sex')
ggsave(paste(outpath,'WindingAPOEProbe2sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp2, Sex %in% "F")
tempMales <- subset(dfp2, Sex %in% "M")

dataF.lm <- lm(Winding ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(Winding ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(Winding ~ APOE, data = tempFemales)
testMethodM<-oneway.test(Winding ~ APOE, data = tempMales)

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'windingAPOEProbe2sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Total Distance in Probe Day 1 for APOE22/33/44
ggline(dfp1, x='APOE', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOEProbe1.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(Distance ~ APOE, data = dfp1)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ APOE, data = dfp1)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp1))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOEProbe1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Total Distance in Probe Day 1 for APOE22/33/44, separated by sex
ggline(dfp1, x='APOE', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolDistAPOEProbe1sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp1, Sex %in% "F")
tempMales <- subset(dfp1, Sex %in% "M")

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOEProbe1sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Total Distance in Probe Day 2 for APOE22/33/44
ggline(dfp2, x='APOE', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOEProbe2.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(Distance ~ APOE, data = dfp2)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ APOE, data = dfp2)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp2))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOEProbe2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Total Distance in Probe Day 2 for APOE22/33/44, separated by sex
ggline(dfp2, x='APOE', y='Distance', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolDistAPOEProbe2sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp2, Sex %in% "F")
tempMales <- subset(dfp2, Sex %in% "M")

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'poolDistAPOEProbe2sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Normalized SW Distance in Probe Day 1 for APOE22/33/44
ggline(dfp1, x='APOE', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top')
ggsave(paste(outpath,'NormSWDistAPOEProbe1.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(NormSWDist ~ APOE, data = dfp1)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(NormSWDist ~ APOE, data = dfp1)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp1))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEProbe1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Normalized SW Distance in Probe Day 1 for APOE22/33/44, separated by sex
ggline(dfp1, x='APOE', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWDistAPOEProbe1sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp1, Sex %in% "F")
tempMales <- subset(dfp1, Sex %in% "M")

dataF.lm <- lm(NormSWDist ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(NormSWDist ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(NormSWDist ~ APOE, data = tempFemales)
testMethodM<-oneway.test(NormSWDist ~ APOE, data = tempMales)

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEProbe1sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Normalized SW Distance in Probe Day 2 for APOE22/33/44
ggline(dfp2, x='APOE', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top')
ggsave(paste(outpath,'NormSWDistAPOEProbe2.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

data.lm <- lm(NormSWDist ~ APOE, data = dfp2)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(NormSWDist ~ APOE, data = dfp2)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfp2))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=5, ncol=5)
postHocTable[1,]=c('', '', '', '', '')
postHocTable[2,]=c('TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTable[3,]=c('APOE3/3-APOE2/2', tukey.test$APOE[1,1], tukey.test$APOE[1,2], tukey.test$APOE[1,3], tukey.test$APOE[1,4])
postHocTable[4,]=c('APOE4/4-APOE2/2', tukey.test$APOE[2,1], tukey.test$APOE[2,2], tukey.test$APOE[2,3], tukey.test$APOE[2,4])
postHocTable[5,]=c('APOE4/4-APOE3/3', tukey.test$APOE[3,1], tukey.test$APOE[3,2], tukey.test$APOE[3,3], tukey.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEProbe2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Normalized SW Distance in Probe Day 2 for APOE22/33/44, separated by sex
ggline(dfp2, x='APOE', y='NormSWDist', color='APOE', fill='APOE',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWDistAPOEProbe2sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempFemales <- subset(dfp2, Sex %in% "F")
tempMales <- subset(dfp2, Sex %in% "M")

dataF.lm <- lm(NormSWDist ~ APOE, data = tempFemales)
dataF.aov <- aov(dataF.lm)
tukeyF.test <- TukeyHSD(dataF.aov)

dataM.lm <- lm(NormSWDist ~ APOE, data = tempMales)
dataM.aov <- aov(dataM.lm)
tukeyM.test <- TukeyHSD(dataM.aov)

testMethodF<-oneway.test(NormSWDist ~ APOE, data = tempFemales)
testMethodM<-oneway.test(NormSWDist ~ APOE, data = tempMales)

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
postHocTableF[3,]=c('APOE3/3-APOE2/2', tukeyF.test$APOE[1,1], tukeyF.test$APOE[1,2], tukeyF.test$APOE[1,3], tukeyF.test$APOE[1,4])
postHocTableF[4,]=c('APOE4/4-APOE2/2', tukeyF.test$APOE[2,1], tukeyF.test$APOE[2,2], tukeyF.test$APOE[2,3], tukeyF.test$APOE[2,4])
postHocTableF[5,]=c('APOE4/4-APOE3/3', tukeyF.test$APOE[3,1], tukeyF.test$APOE[3,2], tukeyF.test$APOE[3,3], tukeyF.test$APOE[3,4])

postHocTableM<-matrix(nrow=5, ncol=5)
postHocTableM[1,]=c('', '', '', '', '')
postHocTableM[2,]=c('Male TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableM[3,]=c('APOE3/3-APOE2/2', tukeyM.test$APOE[1,1], tukeyM.test$APOE[1,2], tukeyM.test$APOE[1,3], tukeyM.test$APOE[1,4])
postHocTableM[4,]=c('APOE4/4-APOE2/2', tukeyM.test$APOE[2,1], tukeyM.test$APOE[2,2], tukeyM.test$APOE[2,3], tukeyM.test$APOE[2,4])
postHocTableM[5,]=c('APOE4/4-APOE3/3', tukeyM.test$APOE[3,1], tukeyM.test$APOE[3,2], tukeyM.test$APOE[3,3], tukeyM.test$APOE[3,4])

myfile<-paste(outpath,'NormSWDistAPOEProbe2sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTableF, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableM, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of Normalized Distance in Quadrant SW in Probe Day 2 for All Genotypes-- PROBE TRIAL
myMean<-aggregate(dfSWp2$DistNorm, by=list(Genotype=dfSWp2$Genotype), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Genotype=dfSWp2$Genotype), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(29)

tabbedMeans<-tapply(myMean$x, list(myMean$Genotype), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Genotype), function(x) c(x=x))
tabbedMeans<-tabbedMeans[7,1:3]
tabbedSE<-tabbedSE[7,1:3]

pdf(file='ProbeNormSWDistP2APOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 2",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Genotype",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file