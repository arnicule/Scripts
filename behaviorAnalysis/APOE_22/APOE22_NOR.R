library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(reshape2)
library(ggpubr)

masterFile='/Users/ar_ni/OneDrive/Desktop/APOE22/NOR/APOE22_NOR.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/APOE22/NOR/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df2<-subset(df, (APOE=='APOE2/2'))
#df3<-subset(df, (APOE=='APOE3/3'))
#df4<-subset(df, (APOE=='APOE4/4'))
dfAll<- rbind(df2)

dfDay1<-subset(dfAll, (Stage=='Dilutions'))
dfDay1_Obj4<-data.frame(dfDay1$Animal, dfDay1$APOE, dfDay1$Sex, dfDay1$Stage, dfDay1$age, dfDay1$Trial, dfDay1$Object4.headtime)
dfDay1_Obj4<-cbind(dfDay1_Obj4, Object='Obj4')
colnames(dfDay1_Obj4)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Age', 'Trial', 'Time', 'Object')

dfDay2<-subset(dfAll, (Stage=='Day2_IdenticalObject' | Stage=='Day2_ShortTermMemory'))
dfDay2_Obj2RI<-data.frame(dfDay2$Animal, dfDay2$APOE, dfDay2$Sex, dfDay2$Stage, dfDay2$age, dfDay2$Trial, dfDay2$Object2_Day2_3_5.RI)
dfDay2_Obj2RI<-cbind(dfDay2_Obj2RI, Object='Obj2')
colnames(dfDay2_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Age', 'Trial', 'RI', 'Object')

dfDay2_Obj2RI_Trial1<-subset(dfDay2_Obj2RI, (Trial=='1'))

dfDay2_Obj2RI_ShortTerm<-subset(dfDay2_Obj2RI, (Trial=='4'))

dfDay3<-subset(dfAll, (Stage=='Day3_Novel Object2'))
dfDay3_Obj2RI<-data.frame(dfDay3$Animal, dfDay3$APOE, dfDay3$Sex, dfDay3$Stage, dfDay3$age, dfDay3$Trial, dfDay3$Object2_Day2_3_5.RI)
dfDay3_Obj2RI<-cbind(dfDay3_Obj2RI, Object='Obj2')
colnames(dfDay3_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Age', 'Trial', 'RI', 'Object')

dfDay5<-subset(dfAll, (Stage=='Day5_NovelObject3'))
dfDay5_Obj2RI<-data.frame(dfDay5$Animal, dfDay5$APOE, dfDay5$Sex, dfDay5$Stage, dfDay5$age, dfDay5$Trial, dfDay5$Object2_Day2_3_5.RI)
dfDay5_Obj2RI<-cbind(dfDay5_Obj2RI, Object='Obj2')
colnames(dfDay5_Obj2RI)<-c('AnimalID', 'Genotype', 'Sex', 'Stage', 'Age', 'Trial', 'RI', 'Object')

dfFin<-rbind(dfDay2_Obj2RI_Trial1, dfDay2_Obj2RI_ShortTerm, dfDay3_Obj2RI, dfDay5_Obj2RI)

#START PLOTTING
#__________________________________________________________________________

#1. Object 4 Head time over Day 1 for APOE 2/2, 3/3, and 4/4 with line
ggplot(dfDay1_Obj4, aes(Trial, Time, fill=Age))+
  
  stat_summary(fun.y=mean, color="black", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_smooth(method='loess', se=FALSE,color='black', aes(group=Age, color=Age))+  #Change se=TRUE to add SE
  ggtitle(paste('Day 1, Object 4 Time'))+
  theme_classic()+
  scale_color_manual(values=c('blueviolet', 'chartreuse1'))+
  scale_fill_manual(values=c('blueviolet', 'chartreuse1'))+
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

ggline(dfDay1_Obj4, x='Trial', y='Time', color='Age', fill='Age',
       error.plot='errorbar', add='mean_se', palette = c('blueviolet', 'chartreuse1'), size=1,
       point.size=1.5, xlab='Trial', ylab='Object 4 Head Time (s)', legend='top')
ggsave(paste(outpath,'Day1_Obj4TimeLineNOR.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#2. Object 4 Head Time over Day 1 for APOE 2/2, 3/3, 4/4
ggerrorplot(dfDay1_Obj4, x='Trial', y='Time', color='Age', fill='Age',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Trial', ylab='Object 4 Head Time (s)',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day1Obj4TimeBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#3. Objects 2 RI over Day 2 
ggerrorplot(dfDay2_Obj2RI, x='Trial', y='RI', color='Age', fill='Age',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Trial', ylab='Object 2 RI',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day2Obj2RIBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)

#4. Object 2 RI over Day 2 LP, Day 2 Short-term Memory, Day 3, and Day 5
ggerrorplot(dfFin, x='Stage', y='RI', color='Age', fill='Age',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='Stage', ylab='Object 2 RI',
            legend='top', position=position_dodge(0.2))
ggsave(paste(outpath, 'Day235Obj2RIBoxNOR.pdf', sep=''), plot=last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)
