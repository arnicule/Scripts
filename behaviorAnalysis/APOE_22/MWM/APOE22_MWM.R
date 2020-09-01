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

dfSWp1<-subset(dfSW, (Day=="Probe_D5"))
dfSWp2<-subset(dfSW, (Day=="Probe_D8"))

dfQuad<-rbind(dfSW, dfSE, dfNW, dfNE)
dfQuad2<-subset(dfQuad, (Genotype=='APOE2/2')) #Info on all quadrants for genotype APOE2/2
dfQuad2<-na.omit(dfQuad2) #Omit n/a
dfQuady<-subset(dfQuad2, (Age=="young"))
dfQuado<-subset(dfQuad2, (Age=="old"))
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

#1. Mean Speed over acquisition day for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='Mean.speed', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Mean Speed (m/s)', legend='top')
ggsave(paste(outpath,'meanSpeedAPOE22.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#1. Mean Speed over acquisition day for APOE 2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='Mean.speed', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Mean Speed (m/s)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'meanSpeedAPOE22Sex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#1. Time to platform over acquisition day for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='Duration', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top')
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

#1. Time to platform over acquisition day for APOE 2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='Duration', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Time to Platform (s)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'poolTimeAPOE22Sex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

temp <- subset(dfAveraged, Stage %in% "Day3")
temp <- subset(temp, Sex %in% "F")

data.lm <- lm(Duration ~ Age, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Duration ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'poolTimeAPOE22Day3Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#2. Distance to Platform over acquisition day for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE22.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

testMethod<-t.test(Distance ~ Age, data = dfAveraged)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Stage %in% "Day5")

testMethod<-t.test(Distance ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfAveraged))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Day5stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
#__________________________________________________________________

#2. Distance to Platform over acquisition day for APOE 2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='', ylab='Distance to Platform (m)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'poolDistAPOE22Sex.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

temp <- subset(dfAveraged, Stage %in% "Day3")
temp <- subset(temp, Sex %in% "F")

data.lm <- lm(Distance ~ Age, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(paste("Female", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Day3Femalestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAveraged, Stage %in% "Day5")
temp <- subset(temp, Sex %in% "M")

data.lm <- lm(Distance ~ Age, data = temp)
data.aov <- aov(data.lm)
tukey.test <- TukeyHSD(data.aov)

testMethod<-oneway.test(Distance ~ Age, data = temp)

mytTable<-as_tibble(
  cbind(paste("Male", testMethod$data.name, sep=" "), testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Day5Malestats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='NormSWDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE22.pdf',sep=''), plot = last_plot(), device = 'pdf',
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

#3. NORMALIZED distance swam in SW quadrant for APOE 2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='NormSWDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Distance (SW)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWDistAPOE22Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, young vs old
ggline(dfAveraged, x='Stage', y='NormSWTime', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top')
ggsave(paste(outpath,'NormSWTimeAPOE22.pdf',sep=''), plot = last_plot(), device = 'pdf',
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

#4. NORMALIZED Time swam in SW quadrant for APOE 2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='NormSWTime', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent of Time (SW)', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWTimeAPOE22Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#5. Normalized SW Time swam adjusted for Mean Speed in APOE2/2, young vs old
ggline(dfAveraged, x='Stage', y='residuals', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top')
ggsave(paste(outpath,'NormSWTimeAdjustAPOE22.pdf',sep=''), plot = last_plot(), device = 'pdf',
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

#5. Normalized SW Time swam adjusted for Mean Speed in APOE2/2, young vs old, separated by sex
ggline(dfAveraged, x='Stage', y='residuals', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent Quadrant Time (SW) - Residuals', legend='top', facet.by="Sex")
ggsave(paste(outpath,'NormSWTimeAdjustAPOE22Sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#_______________________________________________________________________________
#THE REST OF THE SCRIPT MAKES PATTERNED BAR PLOTS FOR PROBE TRIALS
#5. Distance in each Quadrant for young mice in Probes
ggline(dfQuady, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Quadrant Distance (m)', legend='top', facet.by="Day")
ggsave(paste(outpath,'quadDistAPOE22YoungProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Barplot with Standard Error of DISTANCE in Quadrant for age young-- PROBE TRIAL
myMean<-aggregate(dfQuady$Distance, by=list(Day=dfQuady$Day, Quadrant=dfQuady$Quadrant), mean)
mySD<-aggregate(dfQuady$Distance, by=list(Day=dfQuady$Day, Quadrant=dfQuady$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuadyYoungAPOE22.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant APOE2/2 (Young)",
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

dfQuadyp1<-subset(dfQuady, Day=='Probe_D5')
dfQuadyp2<-subset(dfQuady, Day=='Probe_D8')

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

myfile<-paste(outpath,'ProbeQuadDistyoungAPOE22stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#5. Distance in each Quadrant for old mice in Probes
ggline(dfQuado, x='Quadrant', y='Distance', color='Quadrant', fill='Quadrant',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1','red','orange'), size=1, 
       point.size = 1.5, xlab='', ylab='Quadrant Distance (m)', legend='top', facet.by="Day")
ggsave(paste(outpath,'quadDistAPOE22OldProbe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Barplot with Standard Error of DISTANCE in Quadrant for age old-- PROBE TRIAL
myMean<-aggregate(dfQuado$Distance, by=list(Day=dfQuado$Day, Quadrant=dfQuado$Quadrant), mean)
mySD<-aggregate(dfQuado$Distance, by=list(Day=dfQuado$Day, Quadrant=dfQuado$Quadrant), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(11)

tabbedMeans<-tapply(myMean$x, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Quadrant, myMean$Day), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,6:7]
tabbedSE<-tabbedSE[,6:7]

pdf(file='ProbeDistInQuadOldAPOE22.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,.8),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Distance per Quadrant APOE2/2 (Old)",
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

dfQuadop1<-subset(dfQuado, Day=='Probe_D5')
dfQuadop2<-subset(dfQuado, Day=='Probe_D8')

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

myfile<-paste(outpath,'ProbeQuadDistoldAPOE22stats.csv')
write.table(mytTablep1, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTablep2, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTablep1, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTablep2, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#5. SW Distance for all ages in Probes
ggline(dfSW, x='Age', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se',palette = c('blueviolet','chartreuse1'), size=1, 
       point.size = 1.5, xlab='', ylab='SW Quadrant Distance (m)', legend='top', facet.by="Day")
ggsave(paste(outpath,'SWDistAPOE22Probe.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

#Barplot with Standard Error of Distance in Quadrant SW for All ages-- PROBE TRIAL
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

myfile<-paste(outpath,'ProbeQuadDistAPOE22stats.csv')
write.table(mytTabley, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableo, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
write.table(postHocTabley, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)
write.table(postHocTableo, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#Barplot with Standard Error of Normalized Distance in Quadrant SW for All ages in Probe 1-- PROBE TRIAL
myMean<-aggregate(dfSWp1$DistNorm, by=list(Day=dfSWp1$Day, Age=dfSWp1$Age), mean)
mySD<-aggregate(dfSWp1$DistNorm, by=list(Day=dfSWp1$Day, Age=dfSWp1$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(34)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6,1:2]
tabbedSE<-tabbedSE[6,1:2]

pdf(file='ProbeNormSWDistP1APOE22.pdf')
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
                    args.legend=list(title="Age",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW for All ages in Probe 2-- PROBE TRIAL
myMean<-aggregate(dfSWp2$DistNorm, by=list(Day=dfSWp2$Day, Age=dfSWp2$Age), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Day=dfSWp2$Day, Age=dfSWp2$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(34)

tabbedMeans<-tapply(myMean$x, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Day, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[7,1:2]
tabbedSE<-tabbedSE[7,1:2]

pdf(file='ProbeNormSWDistP2APOE22.pdf')
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
                    args.legend=list(title="Age",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE *2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

#Barplot with Standard Error of Normalized Distance in Quadrant SW for All ages in Probe 1-- PROBE TRIAL, separated by sex
myMean<-aggregate(dfSWp1$DistNorm, by=list(Sex=dfSWp1$Sex, Age=dfSWp1$Age), mean)
mySD<-aggregate(dfSWp1$DistNorm, by=list(Sex=dfSWp1$Sex, Age=dfSWp1$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(34)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[6,1:2]
tabbedSE<-tabbedSE[6,1:2]

pdf(file='ProbeNormSWDistP1APOE22Sex.pdf')
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

#Barplot with Standard Error of Normalized Distance in Quadrant SW for All ages in Probe 2-- PROBE TRIAL, separated by sex
myMean<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Age=dfSWp2$Age), mean)
mySD<-aggregate(dfSWp2$DistNorm, by=list(Sex=dfSWp2$Sex, Age=dfSWp2$Age), sd)
myMean<-do.call(data.frame, myMean)
mySD<-do.call(data.frame, mySD)
myMean$SE<-mySD$x/sqrt(34)

tabbedMeans<-tapply(myMean$x, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedSE<-tapply(myMean$SE, list(myMean$Sex, myMean$Age), function(x) c(x=x))
tabbedMeans<-tabbedMeans[7,1:2]
tabbedSE<-tabbedSE[7,1:2]

pdf(file='ProbeNormSWDistP2APOE22Sex.pdf')
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

#Total Distance in Probe Day 1 for APOE22
ggline(dfp1, x='Age', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE22Probe1.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(Distance ~ Age, data = dfp1)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfp1))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Probe1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#Total Distance in Probe Day 1 for APOE22, separated by sex
ggline(dfp1, x='Age', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolDistAPOE22Probe1sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfp1, Sex %in% "F")
tempM <- subset(dfp1, Sex %in% "M")

testMethodF<-t.test(Distance ~ Age, data = tempF)
testMethodM<-t.test(Distance ~ Age, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Probe1Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)

#Total Distance in Probe Day 2 for APOE22
ggline(dfp2, x='Age', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top')
ggsave(paste(outpath,'poolDistAPOE22Probe2.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(Distance ~ Age, data = dfp2)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfp2))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Probe2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#Total Distance in Probe Day 2 for APOE22, separated by sex
ggline(dfp2, x='Age', y='Distance', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Total Distance (m)', legend='top', facet.by='Sex')
ggsave(paste(outpath,'poolDistAPOE22Probe2sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfp2, Sex %in% "F")
tempM <- subset(dfp2, Sex %in% "M")

testMethodF<-t.test(Distance ~ Age, data = tempF)
testMethodM<-t.test(Distance ~ Age, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'poolDistAPOE22Probe2Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)

#Normalized SW Distance in Probe Day 1 for APOE22
ggline(dfp1, x='Age', y='NormDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE22Probe1.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormDist ~ Age, data = dfp1)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfp1))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistAPOE22Probe1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#Normalized SW Distance in Probe Day 1 for APOE22, separated by sex
ggline(dfp1, x='Age', y='NormDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWDistAPOE22Probe1sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfp1, Sex %in% "F")
tempM <- subset(dfp1, Sex %in% "M")

testMethodF<-t.test(NormDist ~ Age, data = tempF)
testMethodM<-t.test(NormDist ~ Age, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistAPOE22Probe1Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)

#Normalized SW Distance in Probe Day 2 for APOE22
ggline(dfp2, x='Age', y='NormDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top')
ggsave(paste(outpath,'NormSWDistAPOE22Probe2.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

testMethod<-t.test(NormDist ~ Age, data = dfp2)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$estimate[1], testMethod$estimate[2], testMethod$statistic, testMethod$p.value,  
        testMethod$conf.int[1], testMethod$conf.int[2] , testMethod$parameter, nrow(dfp2))
)

mycolnames<-c('contrast', names(testMethod$estimate)[1], names(testMethod$estimate)[2], names(testMethod$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDisttAPOE22Probe2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

#Normalized SW Distance in Probe Day 2 for APOE22, separated by sex
ggline(dfp2, x='Age', y='NormDist', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet','chartreuse1','red'), size=1, 
       point.size = 1.5, xlab='', ylab='Percent SW Distance', legend='top', facet.by='Sex')
ggsave(paste(outpath,'NormSWDistAPOE22Probe2sex.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 5, height = 5, units = c("in"),dpi = 300)

tempF <- subset(dfp2, Sex %in% "F")
tempM <- subset(dfp2, Sex %in% "M")

testMethodF<-t.test(NormDist ~ Age, data = tempF)
testMethodM<-t.test(NormDist ~ Age, data = tempM)

mytTableF<-as_tibble(
  cbind(paste("Female", testMethodF$data.name, sep=" "), testMethodF$estimate[1], testMethodF$estimate[2], testMethodF$statistic, testMethodF$p.value,  
        testMethodF$conf.int[1], testMethodF$conf.int[2] , testMethodF$parameter, nrow(tempF))
)
mytTableM<-as_tibble(
  cbind(paste("Male", testMethodM$data.name, sep=" "), testMethodM$estimate[1], testMethodM$estimate[2], testMethodM$statistic, testMethodM$p.value,  
        testMethodM$conf.int[1], testMethodM$conf.int[2] , testMethodM$parameter, nrow(tempM))
)

mycolnames<-c('contrast', names(testMethodF$estimate)[1], names(testMethodF$estimate)[2], names(testMethodF$statistic), 'pvalue', 'CIlwr', 'CIhi', 'df', 'observations')

myfile<-paste(outpath,'NormSWDistAPOE22Probe2Sexstats.csv')
write.table(mytTableF, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(mytTableM, file=myfile, col.names = FALSE, sep = "," , row.names = F,append=TRUE)
