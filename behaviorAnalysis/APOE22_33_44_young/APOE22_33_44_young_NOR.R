library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(reshape2)
library(ggpubr)

masterFile='/Users/ar_ni/OneDrive/Desktop/APOE_NOR/APOE_22_33_44_NORcombined.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/APOE_NOR/R_Graphs/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df2<-subset(df, (APOE=='APOE2/2'))
df3<-subset(df, (APOE=='APOE3/3'))
df4<-subset(df, (APOE=='APOE4/4'))
dfAll<- rbind(df2, df3, df4)

dfOpen<-subset(dfAll, (Stage=='openfield'))
dfOpen<-droplevels(dfOpen)

dfLP<-subset(dfAll, (Stage=='Dilutions') & (Trial=="1"))

dfDay1<-subset(dfAll, (Stage=='Dilutions'))
dfDay1_Obj4<-data.frame(dfDay1$Animal, dfDay1$APOE, dfDay1$Sex, dfDay1$Stage, dfDay1$Trial, dfDay1$Object4.headtime)
dfDay1_Obj4<-cbind(dfDay1_Obj4, Object='Obj4')
colnames(dfDay1_Obj4)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Trial', 'Time', 'Object')

dfDay2<-subset(dfAll, (Stage=='Day2_IdenticalObject' | Stage=='Day2_ShortTermMemory'))
dfDay2_Obj2RI<-data.frame(dfDay2$Animal, dfDay2$APOE, dfDay2$Sex, dfDay2$Stage, dfDay2$Trial, dfDay2$Object2_Day2_3_5.RI)
dfDay2_Obj2RI<-cbind(dfDay2_Obj2RI, Object='Obj2')
colnames(dfDay2_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Trial', 'RI', 'Object')

dfDay2_Obj2RI_ShortTerm<-subset(dfDay2_Obj2RI, (Trial=='4'))

dfDay3<-subset(dfAll, (Stage=='Day3_Novel Object2'))
dfDay3_Obj2RI<-data.frame(dfDay3$Animal, dfDay3$APOE, dfDay3$Sex, dfDay3$Stage, dfDay3$Trial, dfDay3$Object2_Day2_3_5.RI)
dfDay3_Obj2RI<-cbind(dfDay3_Obj2RI, Object='Obj2')
colnames(dfDay3_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Trial', 'RI', 'Object')

dfDay5<-subset(dfAll, (Stage=='Day5_NovelObject3'))
dfDay5_Obj2RI<-data.frame(dfDay5$Animal, dfDay5$APOE, dfDay5$Sex, dfDay5$Stage, dfDay5$Trial, dfDay5$Object2_Day2_3_5.RI)
dfDay5_Obj2RI<-cbind(dfDay5_Obj2RI, Object='Obj2')
colnames(dfDay5_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Trial', 'RI', 'Object')

dfFin<-rbind(dfDay2_Obj2RI_ShortTerm, dfDay3_Obj2RI, dfDay5_Obj2RI)

#START PLOTTING
#__________________________________________________________________________

#1. Object 4 Head time over Day 1 for APOE 2/2, 3/3, and 4/4 with line
ggplot(dfDay1_Obj4, aes(Trial, Time, fill=Genotype))+
  
  stat_summary(fun.y=mean, color="black", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_smooth(method='loess', se=FALSE,color='black', aes(group=Genotype, color=Genotype))+  #Change se=TRUE to add SE
  ggtitle(paste('Day 1, Object 4 Time'))+
  theme_classic()+
  scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  scale_fill_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  theme(legend.position="top")+
  labs(y='Time (sec)', x='Trial')+
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())

ggline(dfDay1_Obj4, x='Trial', y='Time', color='Genotype', fill='Genotype',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='Trial', ylab='Object 4 Head Time (s)', legend='top')
ggsave(paste(outpath,'Day1_Obj4TimeLineNOR.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#2. Object 4 Head Time over Day 1 for APOE 2/2, 3/3, 4/4
ggerrorplot(dfDay1_Obj4, x='Trial', y='Time', color='Genotype', fill='Genotype',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1', 'red'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Trial', ylab='Object 4 Head Time (s)',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day1Obj4TimeBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#3. Objects 2 RI over Day 2 
ggerrorplot(dfDay2_Obj2RI, x='Trial', y='RI', color='Genotype', fill='Genotype',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1', 'red'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Trial', ylab='Object 2 RI',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day2Obj2RIBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

temp <- subset(dfDay2_Obj2RI, Trial %in% "1")

testMethod<-oneway.test(RI ~ Genotype, data = temp)
postHoc<-pairwise.t.test(temp$RI, temp$Genotype)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=4, ncol=3)
postHocTable[1,]=c('', '', '')
postHocTable[2,]=c('Pairwise Comparisons (p-value)', 'APOE2/2', 'APOE3/3')
postHocTable[3,]=c('APOE3/3', postHoc$p.value[1,1], postHoc$p.value[1,2])
postHocTable[4,]=c('APOE4/4', postHoc$p.value[2,1], postHoc$p.value[2,2])

myfile<-paste(outpath,'Obj2Day2Trial1stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfDay2_Obj2RI, Trial %in% "2")

testMethod<-oneway.test(RI ~ Genotype, data = temp)
postHoc<-pairwise.t.test(temp$RI, temp$Genotype)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=4, ncol=3)
postHocTable[1,]=c('', '', '')
postHocTable[2,]=c('Pairwise Comparisons (p-value)', 'APOE2/2', 'APOE3/3')
postHocTable[3,]=c('APOE3/3', postHoc$p.value[1,1], postHoc$p.value[1,2])
postHocTable[4,]=c('APOE4/4', postHoc$p.value[2,1], postHoc$p.value[2,2])

myfile<-paste(outpath,'Obj2Day2Trial2stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

temp <- subset(dfDay2_Obj2RI, Trial %in% "4")

testMethod<-oneway.test(RI ~ Genotype, data = temp)
postHoc<-pairwise.t.test(temp$RI, temp$Genotype)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

postHocTable<-matrix(nrow=4, ncol=3)
postHocTable[1,]=c('', '', '')
postHocTable[2,]=c('Pairwise Comparisons (p-value)', 'APOE2/2', 'APOE3/3')
postHocTable[3,]=c('APOE3/3', postHoc$p.value[1,1], postHoc$p.value[1,2])
postHocTable[4,]=c('APOE4/4', postHoc$p.value[2,1], postHoc$p.value[2,2])

myfile<-paste(outpath,'Obj2Day2Trial4stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(postHocTable, file=myfile, sep=",", row.names=F, append=TRUE, col.names=F)

#4. Object 2 RI over Day 2 Short-term Memory, Day 3, and Day 5
ggerrorplot(dfFin, x='Stage', y='RI', color='Genotype', fill='Genotype',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1', 'red'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Stage', ylab='Object 2 RI',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day235Obj2RIBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)



#Location Preferences
#9. Barplot with Standard Error of Head times for all objects in Day 1 Trial 1
myLP<-subset(dfLP, select=c(APOE, Object1.headtime, Object2.headtime, Object.3.headtime, Object4.headtime))

myMean<-aggregate(myLP, by=list(APOE=myLP$APOE), mean)
myMean<-myMean[-c(2)]
myMean.long<-melt(myMean,id.vars="APOE")
mySD<-aggregate(myLP, by=list(APOE=myLP$APOE), sd)
mySD<-mySD[-c(2)]
mySD.long<-melt(mySD,id.vars="APOE")
myMean.long<-do.call(data.frame, myMean.long)
mySD.long<-do.call(data.frame, mySD.long)
myMean.long$SE<-mySD.long$value/sqrt(9)

head(subset(myLP, (myLP$APOE==mySD.long$APOE)), 25)
head(myMean.long, 25)
head(mySD.long, 25)

tabbedMeans<-tapply(myMean.long$value, list(myMean.long$variable, myMean.long$APOE), function(x) c(x=x))
tabbedSE<-tapply(myMean.long$SE, list(myMean.long$variable, myMean.long$APOE), function(x) c(x=x))
tabbedMeans<-tabbedMeans[,1:3]
tabbedSE<-tabbedSE[,1:3]

pdf(file='LPAPOE.pdf')
barCenters<-barplot(height=tabbedMeans,
                    density=c(100,5,15,30),
                    angle=c(0,0,45,90),
                    ylim=c(0,14),
                    col='blueviolet',
                    beside=TRUE, las=1,
                    cex.names=0.75,
                    main="Head time",
                    ylab="Time (s)",
                    border="black", axes=TRUE,
                    legend.text=TRUE,
                    args.legend=list(title="Object",
                                     x="topright",
                                     cex=.7))
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE * 2, lwd = 1.5)
arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off() #Close pdf file

# Boxplots for Centerzone and Outerzone head time of all genotypes
boxplot(Centerzone.headtime ~ APOE, data = dfOpen, ylim=range(0:300))
ggsave(paste(outpath, 'openfieldcentertimeAPOE.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

boxplot(Outerzone.headtime ~ APOE, data = dfOpen, ylim=range(100:500))
ggsave(paste(outpath, 'openfieldoutertimeAPOE.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

# Boxplots for Centerzone and Outerzone distance of all genotypes
boxplot(Centerzone.distance ~ APOE, data = dfOpen, ylim=range(0:10))
ggsave(paste(outpath, 'openfieldcenterdistAPOE.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

boxplot(Outerzone.distance ~ APOE, data = dfOpen, ylim=range(0:30))
ggsave(paste(outpath, 'openfieldouterdistAPOE.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


#9. Barplot with Standard Error of Head times for all objects in Day 1 Trial 1
ggline(dfDay1_Obj4, x='Trial', y='Time', color='Genotype', fill='Genotype',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1', 'red'), size=1,
       point.size=1.5, xlab='Trial', ylab='Object 4 Head Time (s)', legend='top')