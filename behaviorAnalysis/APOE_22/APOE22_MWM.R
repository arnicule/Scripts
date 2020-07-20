library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(ggpubr)

masterFile='/Users/ar_ni/OneDrive/Desktop/APOE22/MWM/MWMDataAPOE22.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/APOE22/MWM/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df2<-subset(df, (APOE=='APOE2/2')) #Keeps only Genotypes 2/2, 3/3, and 4/4
dfAll<-rbind(df2)

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
dfAveraged<-aggregate(.~Animal+APOE+Sex+Age+Stage, dfFin, mean, na.action=na.pass)
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
dfNW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Age, dfProbe$NW.time, dfProbe$NW.distance, dfProbe$NW.Time.Norm, dfProbe$NW.Dist.Norm)
dfNW<-cbind(dfNW, quadrant='NW')
colnames(dfNW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Age', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSW<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Age, dfProbe$SW.time, dfProbe$SW.distance, dfProbe$SW.Time.Norm, dfProbe$SW.Dist.Norm)
dfSW<-cbind(dfSW, quadrant='SW')
colnames(dfSW)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Age', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfNE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Age, dfProbe$NE.time, dfProbe$NE.distance, dfProbe$NE.Time.Norm, dfProbe$NE.Dist.Norm)
dfNE<-cbind(dfNE, quadrant='NE')
colnames(dfNE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Age', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')
dfSE<-data.frame(dfProbe$Animal, dfProbe$APOE, dfProbe$Stage, dfProbe$Sex, dfProbe$Age, dfProbe$SE.time, dfProbe$SE.distance, dfProbe$SE.Time.Norm, dfProbe$SE.Dist.Norm)
dfSE<-cbind(dfSE, quadrant='SE')
colnames(dfSE)<-c('Animal', 'Genotype', 'Day', 'Sex', 'Age', 'Time', 'Distance', 'TimeNorm', 'DistNorm', 'Quadrant')

dfQuad<-rbind(dfSW, dfSE, dfNW, dfNE)
dfQuad2<-subset(dfQuad, (Genotype=='APOE2/2')) #Info on all quadrants for genotype APOE2/2
#dfQuad3<-subset(dfQuad, (Genotype=='APOE3/3')) #Info on all quadrants for genotype APOE3/3
#dfQuad4<-subset(dfQuad, (Genotype=='APOE4/4')) #Info on all quadrants for genotype APOE4/4
dfQuad2<-na.omit(dfQuad2) #Omit n/a
#dfQuad3<-na.omit(dfQuad3) #Omit n/a
#dfQuad4<-na.omit(dfQuad4) #Omit n/a
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

#1. Time to platform over acquisition day for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='Time', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Mean Speed', legend='top')
ggsave(paste(outpath,'poolTimeAPOE22.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

testMethod<-t.test(Duration ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolTimeAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________________________

#2. Distance to Platform over acquisition day for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

testMethod<-t.test(Distance ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________________________

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='NormSWDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormSWDist ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Stage %in% "Day1")

testMethod<-t.test(NormSWDist ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistAPOE22Day1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#___________________________________________________________________
#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='NormSWTime', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top')
ggsave(paste(outpath,'NormSWTimeAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormSWTime ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWTimeAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#_________________________________________________________________________
#5. SW Time swam adjusted for Mean Speed in APOE2/2, young vs old
ggline(dfAveraged, x='Stage', y='residuals', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'SWTimeAdjustAPOE.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(residuals ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'SWTimeAdjustAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <-subset(dfAveraged, Stage %in% "Day1")

testMethod<-t.test(residuals ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'SWTimeAdjustAPOE22Day1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <-subset(dfAveraged, Stage %in% "Day2")

testMethod<-t.test(residuals ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'SWTimeAdjustAPOE22Day2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS
#6. Barplot with Standard Error of DISTANCE in Quadrant for Genotype 2/2-- PROBE TRIAL
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
myMean<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Age=dfSW$Age), mean)
mySD<-aggregate(dfSW$Distance, by=list(Day=dfSW$Day, Age=dfSW$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6:7,1:2]
tabbedSE<-tabbedSE[6:7,1:2]

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