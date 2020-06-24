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
dfNW<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$Day, dfProbe$NW.Time, dfProbe$NW.Distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW,quadrant='NW')
colnames(dfNW)<-c('Animal', 'genotype', 'Age','Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$Day, dfProbe$SW.Time, dfProbe$SW.Distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW,quadrant='SW')
colnames(dfSW)<-c('Animal', 'genotype', 'Age','Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$Day, dfProbe$NE.Time, dfProbe$NE.Distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE,quadrant='NE')
colnames(dfNE)<-c('Animal', 'genotype', 'Age','Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$genotype, dfProbe$age_group, dfProbe$Day, dfProbe$SE.Time, dfProbe$SE.Distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE,quadrant='SE')
colnames(dfSE)<-c('Animal', 'genotype', 'Age','Day', 'Time','Distance','TimeNorm', 'DistNorm','Quadrant')

dfQuad<-rbind(dfSW,dfSE,dfNW,dfNE)
dfQuady<-subset(dfQuad, (Age=='young'))
dfQuado<-subset(dfQuad, (Age=='old'))
dfQuady<-na.omit(dfQuady)
dfQuado<-na.omit(dfQuado)
dfQuadp1<-subset(dfQuad, Day=='ProbeTrial1')
dfQuadp2<-subset(dfQuad, Day=='ProbeTrial2')

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
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'poolTimeC57stats.csv')
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
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'poolDistC57stats.csv')
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

#________________________________________________________________________



#4. Time spend in SW Quadrant for ages young vs. old separated by sex
ggline(dfAveraged, x='Day', y='SW.Time', color='age_group', fill='age_group',
          error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
          point.size=1.5, xlab='', ylab='SW Time (s)', legend='top',
          facet.by="sex"
)
ggsave(paste(outpath,'SWTimeC57Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

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
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________
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
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

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
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(dfAveraged))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'MeanSpeedC57stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
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