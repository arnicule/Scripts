library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(ggpubr)

masterFile='/Users/ar_ni/OneDrive/Desktop/C57_MWM/mwm_master_organized_MWM.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/C57_MWM/R_Graphs/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df0<-subset(df, (genotype=='0')) #Keeps only genotype 0 (C57)
dfAll<-df0
dfAll$age_group<-as.factor(dfAll$age_group) #Converts age_group number to factor
head(dfAll)

#Correct for inaccurate distance calibration in AnyMaze
dfAll$Pool.Distance=dfAll$Pool.Distance/10
dfAll$NE.Distance=dfAll$NE.Distance/10
dfAll$NW.Distance=dfAll$NW.Distance/10
dfAll$SE.Distance=dfAll$SE.Distance/10
dfAll$SW.Distance=dfAll$SW.Distance/10

#Remove Probe Trials from Data Set
df1<-subset(dfAll, (Day=='Day1'))
df2<-subset(dfAll, (Day=='Day2'))
df3<-subset(dfAll, (Day=='Day3'))
df4<-subset(dfAll, (Day=='Day4'))
df5<-subset(dfAll, (Day=='Day5'))
dfFin<-rbind(df1,df2,df3,df4,df5) #dfFin contains info from all regular trials

#Normalize time and distance in target region
dfFin$NormSWTime<-dfFin$SW.Time/dfFin$Pool.time
dfFin$NormSWDist<-dfFin$SW.Distance/dfFin$Pool.Distance
dfFin<-subset(dfFin, (NormSWTime <= 1))

#Averages 4 trials per day for each mouse
dfAveraged<-aggregate(.~Animal+genotype+sex+runno+Day+age_group, dfFin, mean, na.action=na.pass)

#Extract probe trials only 
dfp1<-subset(dfAll, (Day=='ProbeTrial1'))
dfp2<-subset(dfAll, (Day=='ProbeTrial2'))
dfProbe<-rbind(dfp1,dfp2) #dfProbe contains info from all probe trials

#Normalize Probe Distances
dfProbe$DistTot<-dfProbe$NE.Distance+dfProbe$NW.Distance+dfProbe$SE.Distance+dfProbe$SW.Distance
dfProbe$NE.Dist.Norm<-dfProbe$NE.Distance/dfProbe$DistTot
dfProbe$NW.Dist.Norm<-dfProbe$NW.Distance/dfProbe$DistTot
dfProbe$SE.Dist.Norm<-dfProbe$SE.Distance/dfProbe$DistTot
dfProbe$SW.Dist.Norm<-dfProbe$SW.Distance/dfProbe$DistTot

#Normalize Probe Times
dfProbe$TimeTot<-dfProbe$NE.Time+dfProbe$NW.Time+dfProbe$SE.Time+dfProbe$SW.Time
dfProbe$NE.Time.Norm<-dfProbe$NE.Time/dfProbe$TimeTot
dfProbe$NW.Time.Norm<-dfProbe$NW.Time/dfProbe$TimeTot
dfProbe$SE.Time.Norm<-dfProbe$SE.Time/dfProbe$TimeTot
dfProbe$SW.Time.Norm<-dfProbe$SW.Time/dfProbe$TimeTot

#This section adds separates information by quadrants and adds quadrant labels
dfNW<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$sex, dfProbe$Day, dfProbe$NW.Time, dfProbe$NW.Distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW,quadrant='NW')
colnames(dfNW)<-c('Animal', 'genotype', 'Age', 'Sex', 'Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$sex, dfProbe$Day, dfProbe$SW.Time, dfProbe$SW.Distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW,quadrant='SW')
colnames(dfSW)<-c('Animal', 'genotype', 'Age', 'Sex', 'Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$sex, dfProbe$Day, dfProbe$NE.Time, dfProbe$NE.Distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE,quadrant='NE')
colnames(dfNE)<-c('Animal', 'genotype', 'Age', 'Sex', 'Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$sex, dfProbe$Day, dfProbe$SE.Time, dfProbe$SE.Distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE,quadrant='SE')
colnames(dfSE)<-c('Animal', 'genotype', 'Age', 'Sex', 'Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')

dfSWp1<-subset(dfSW, (Day=='ProbeTrial1'))
dfSWp2<-subset(dfSW, (Day=='ProbeTrial2'))

dfQuad<-rbind(dfSW,dfSE,dfNW,dfNE)
dfQuady<-subset(dfQuad, (Age=='young'))
dfQuado<-subset(dfQuad, (Age=='old'))
dfQuady<-na.omit(dfQuady)
dfQuado<-na.omit(dfQuado)
dfQuadp1<-subset(dfQuad, Day=='ProbeTrial1')
dfQuadp2<-subset(dfQuad, Day=='ProbeTrial2')

#Adjust SW Time for covariate Mean Speed
dfAveragedtemp<-dfAveraged[-c(22,23,42,61,63,121,122),]
cor.test(dfAveragedtemp$NormSWTime, dfAveragedtemp$Average.Pool.Speed)
lm1<-lm(NormSWTime ~ Average.Pool.Speed, data=dfAveragedtemp) #Na.action not working!
cor.test(lm1$residuals, dfAveragedtemp$Average.Pool.Speed)
dfAveragedtemp<-cbind(dfAveragedtemp, residuals=lm1$residuals)

#__________________________________________________________________________
#Rewrite Data frame for Stats
write.csv(dfAveraged, file='mwmStatsAvgC57.csv')
write.csv(dfQuadp1, file='mwmProbe1C57.csv')
write.csv(dfQuadp2, file='mwmProbe2C57.csv')


#START PLOTTING
#________________________________________________________________________

#1. Time to platform for ages young vs. old
ggline(dfAveraged, x='Day', y='Pool.time', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top',
       )
ggsave(paste(outpath,'poolTimeC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(Pool.time ~ age_group, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Day %in% "Day1")

testMethod<-t.test(Pool.time ~ age_group, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeC57Day1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#________________________________________________________________________

#2. Distance to platform for ages young vs. old
ggline(dfAveraged, x='Day', y='Pool.Distance', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top',
)
ggsave(paste(outpath,'poolDistC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(Pool.Distance ~ age_group, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Day %in% "Day1")

testMethod<-t.test(Pool.Distance ~ age_group, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistC57Day1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#________________________________________________________________________

#3. Distance to platform for ages young vs. old separated by sex
ggline(dfAveraged, x='Day', y='Pool.Distance', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top',
       facet.by='sex'
)
ggsave(paste(outpath,'poolDistC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfAveraged, sex %in% "F")
tempM <- subset(dfAveraged, sex %in% "M")

testMethodF<-t.test(Pool.Distance ~ age_group, data = tempF)
testMethodM<-t.test(Pool.Distance ~ age_group, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistC57Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Day %in% "Day1")
temp <- subset(temp, sex %in% "F")

testMethod<-t.test(Pool.Distance ~ age_group, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistC57Day1Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#________________________________________________________________________
#4. Time spend in Pool for ages young vs. old separated by sex
ggline(dfAveraged, x='Day', y='Pool.time', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top',
       facet.by="sex"
)
ggsave(paste(outpath,'PoolTimeC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfAveraged, sex %in% "F")
tempM <- subset(dfAveraged, sex %in% "M")

testMethodF<-t.test(Pool.time ~ age_group, data = tempF)
testMethodM<-t.test(Pool.time ~ age_group, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeC57Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, sex %in% "F")
temp <- subset(dfAveraged, Day %in% "Day1")

testMethod<-t.test(Pool.time ~ age_group, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(temp))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeC57Day1Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#________________________________________________________________________
#4. Normalized time spend in SW Quadrant for ages young vs. old separated by sex
ggline(dfAveraged, x='Day', y='NormSWTime', color='age_group', fill='age_group',
          error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
          point.size=1.5, xlab='', ylab='Percent SW Time', legend='top',
          facet.by="sex"
)
ggsave(paste(outpath,'NormSWTimeC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfAveraged, sex %in% "F")
tempM <- subset(dfAveraged, sex %in% "M")

testMethodF<-t.test(SW.Time ~ age_group, data = tempF)
testMethodM<-t.test(SW.Time ~ age_group, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWTimeC57Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
#________________________________________________________________________

#5. NORMALIZED distance swam in SW quadrant for ages young vs old
ggline(dfAveraged, x='Day', y='NormSWDist', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance', legend='top', title='Normalized SW Distance (C57)',
       #facet.by="sex" #Uncomment to group x-axis by sex
       )
ggsave(paste(outpath,'NormSWDistC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormSWDist ~ age_group, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________

#5. NORMALIZED distance swam in SW quadrant for ages young vs old, separated by sex
ggline(dfAveraged, x='Day', y='NormSWDist', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Percent of Distance', legend='top', title='Normalized SW Distance (C57)',
       facet.by="sex"
)
ggsave(paste(outpath,'NormSWDistC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#6. NORMALIZED Time swam in SW quadrant for ages young vs old
ggline(dfAveraged, x='Day', y='NormSWTime', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time', legend='top', title='Normalized SW Time (C57)',
       #facet.by="sex" #Uncomment to group x-axis by sex
       )
ggsave(paste(outpath,'NormSWTimeC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormSWTime ~ age_group, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWTimeC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#_______________________________________________________________________________

#7. Mean speed for ages young vs old
ggline(dfAveraged, x='Day', y='Average.Pool.Speed', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', title='Mean Speed (C57)',
)
ggsave(paste(outpath,'MeanSpeedC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(Average.Pool.Speed ~ age_group, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'MeanSpeedC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Day %in% "Day1")

testMethod<-t.test(Average.Pool.Speed ~ age_group, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'MeanSpeedC57Day1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#7. Mean speed for ages young vs old, separated by sex
ggline(dfAveraged, x='Day', y='Average.Pool.Speed', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', title='Mean Speed (C57)', facet.by="sex"
)
ggsave(paste(outpath,'MeanSpeedC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#7. Normalized SW Time, adjusted for speed, for ages young vs old
ggline(dfAveragedtemp, x='Day', y='residuals', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Time - Residuals', legend='top', title='NormSWTime Adjusted (C57)'
)
ggsave(paste(outpath,'NormSWTimeAdjustC57.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#7. Normalized SW Time, adjusted for speed, for ages young vs old, separated by sex
ggline(dfAveragedtemp, x='Day', y='residuals', color='age_group', fill='age_group',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet', 'chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Time - Residuals', legend='top', title='NormSWTime Adjusted (C57)', facet.by='sex'
)
ggsave(paste(outpath,'NormSWTimeAdjustC57sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS
#8. Barplot with Standard Error of DISTANCE in Quadrant for Young mice-- Probe Trial
myMean<-aggregate(dfQuady$Distance, by=list(Day=dfQuady$Day, Quadrant=dfQuady$Quadrant), mean)
mySD<-aggregate(dfQuady$Distance, by=list(Day=dfQuady$Day, Quadrant=dfQuady$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuad_youngC57.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,6),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant (Young)",
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

dfQuadyp1<-subset(dfQuady, Day=='ProbeTrial1')
dfQuadyp2<-subset(dfQuady, Day=='ProbeTrial2')

datap1.lm <- lm(Distance ~ Quadrant, data = dfQuadyp1)
datap1.aov <- aov(datap1.lm)
tukeyp1.test <- TukeyHSD(datap1.aov)
tukeyp1.test

datap2.lm <- lm(Distance ~ Quadrant, data = dfQuadyp2)
datap2.aov <- aov(datap2.lm)
tukeyp2.test <- TukeyHSD(datap2.aov)

testMethodp1<-oneway.test(Distance ~ Quadrant, data = dfQuadyp1)
testMethodp2<-oneway.test(Distance ~ Quadrant, data = dfQuadyp2)

mytTablep1<-as_tibble(
  cbind(paste("Probe 1", testMethodp1$data.name, sep=" "), testMethodp1$statistic, testMethodp1$p.value, testMethodp1$parameter[1], nrow(dfQuadyp1)) #Get values from summary
)

mytTablep2<-as_tibble(
  cbind(paste("Probe 2", testMethodp2$data.name, sep=" "), testMethodp2$statistic, testMethodp2$p.value, testMethodp2$parameter[1], nrow(dfQuadyp2)) #Get values from summary
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

myfile<-paste(outpath,'ProbeQuadDistyoungC57stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of DISTANCE in Quadrant for Old mice-- Probe Trial 
myMean<-aggregate(dfQuado$Distance, by=list(Day=dfQuado$Day, Quadrant=dfQuado$Quadrant), mean)
mySD<-aggregate(dfQuado$Distance, by=list(Day=dfQuado$Day, Quadrant=dfQuado$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuad_oldC57.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,6),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant (Old)",
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

dfQuadop1<-subset(dfQuado, Day=='ProbeTrial1')
dfQuadop2<-subset(dfQuado, Day=='ProbeTrial2')

datap1.lm <- lm(Distance ~ Quadrant, data = dfQuadop1)
datap1.aov <- aov(datap1.lm)
tukeyp1.test <- TukeyHSD(datap1.aov)
tukeyp1.test

datap2.lm <- lm(Distance ~ Quadrant, data = dfQuadop2)
datap2.aov <- aov(datap2.lm)
tukeyp2.test <- TukeyHSD(datap2.aov)

testMethodp1<-oneway.test(Distance ~ Quadrant, data = dfQuadop1)
testMethodp2<-oneway.test(Distance ~ Quadrant, data = dfQuadop2)

mytTablep1<-as_tibble(
  cbind(paste("Probe 1", testMethodp1$data.name, sep=" "), testMethodp1$statistic, testMethodp1$p.value, testMethodp1$parameter[1], nrow(dfQuadop1)) #Get values from summary
)

mytTablep2<-as_tibble(
  cbind(paste("Probe 2", testMethodp2$data.name, sep=" "), testMethodp2$statistic, testMethodp2$p.value, testMethodp2$parameter[1], nrow(dfQuadop2)) #Get values from summary
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

myfile<-paste(outpath,'ProbeQuadDistoldC57stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of Distance in Quadrant SW for both ages-- PROBE TRIAL
myMean<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Age=dfSW$Age), mean)
mySD<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Age=dfSW$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6:7,1:2]
tabbedSE<-tabbedSE[6:7,1:2]

pdf(file='ProbeDistInQuadSWC57.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,6),
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

dfQuadySW<-subset(dfSW, Age=='young')
dfQuadoSW<-subset(dfSW, Age=='old')

datay.lm <- lm(Distance ~ Day, data = dfQuadySW)
datay.aov <- aov(datay.lm)
tukeyy.test <- TukeyHSD(datay.aov)

datao.lm <- lm(Distance ~ Day, data = dfQuadoSW)
datao.aov <- aov(datao.lm)
tukeyo.test <- TukeyHSD(datao.aov)

testMethody<-oneway.test(Distance ~ Day, data = dfQuadySW)
testMethodo<-oneway.test(Distance ~ Day, data = dfQuadoSW)

mytTabley<-as_tibble(
  cbind(paste("Young", testMethody$data.name, sep=" "), testMethody$statistic, testMethody$p.value, testMethody$parameter[1], nrow(dfQuadySW)) #Get values from summary
)
mytTableo<-as_tibble(
  cbind(paste("Old", testMethodo$data.name, sep=" "), testMethodo$statistic, testMethodo$p.value, testMethodo$parameter[1], nrow(dfQuadoSW)) #Get values from summary
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTabley<-matrix(nrow=3, ncol=5)
postHocTabley[1,]=c('', '', '', '', '')
postHocTabley[2,]=c('Young TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTabley[3,]=c('Probe2-Probe1', tukeyy.test$Day[1,1], tukeyy.test$Day[1,2], tukeyy.test$Day[1,3], tukeyy.test$Day[1,4])

postHocTableo<-matrix(nrow=3, ncol=5)
postHocTableo[1,]=c('', '', '', '', '')
postHocTableo[2,]=c('Old TukeyHSD', 'mean diff', 'CIlwr', 'CIhi', 'p-value')
postHocTableo[3,]=c('Probe2-Probe1', tukeyo.test$Day[1,1], tukeyo.test$Day[1,2], tukeyo.test$Day[1,3], tukeyo.test$Day[1,4])

myfile<-paste(outpath,'ProbeQuadDistC57stats.csv')
write.table(mytTabley, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableo, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTabley, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableo, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of Normalized Distance in Quadrant SW for both ages in Probe 1-- PROBE TRIAL
dfSWp1temp<-dfSWp1[-c(16,23),]   #Remove rows with NaN
myMean<-aggregate(dfSWp1temp$DistNorm, by=list(Age=dfSWp1temp$Age), mean)
mySD<-aggregate(dfSWp1temp$DistNorm, by=list(Age=dfSWp1temp$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(28)

tabbedMeans<-tapply(myMean$x, list(myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Age), function(x) c(x=x))

pdf(file='ProbeNormSWDistP1C57.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,0.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 1",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Age",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW for both ages in Probe 2-- PROBE TRIAL
myMean<-aggregate(dfSWp2$DistNorm, by=list(Age=dfSWp2$Age), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Age=dfSWp2$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(30)

tabbedMeans<-tapply(myMean$x, list(myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[7,1:2]
tabbedSE<-tabbedSE[7,1:2]

pdf(file='ProbeNormSWDistP2C57.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,0.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Percent Distance in Quadrant SW Probe 2",
                    ylab="Percent Distance in SW",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Age",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW for both ages in Probe 1-- PROBE TRIAL
#separated by sex
dfSWp1temp<-dfSWp1[-c(16,23),]   #Remove rows with NaN
myMean<-aggregate(dfSWp1temp$DistNorm, by=list(Sex=dfSWp1temp$Sex, Age=dfSWp1temp$Age), mean)
mySD<-aggregate(dfSWp1temp$DistNorm, by=list(Sex=dfSWp1temp$Sex, Age=dfSWp1temp$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(28)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Age), function(x) c(x=x))

pdf(file='ProbeNormSWDistP1C57Sex.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,0.8),
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

#Barplot with Standard Error of Normalized Distance in Quadrant SW for both ages in Probe 2-- PROBE TRIAL
#separated by sex
myMean<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Age=dfSWp2$Age), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Age=dfSWp2$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(30)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Age), function(x) c(x=x))

pdf(file='ProbeNormSWDistP2C57Sex.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,0.8),
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
