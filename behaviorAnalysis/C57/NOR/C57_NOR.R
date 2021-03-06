library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(reshape2)

masterFile='/Users/ar_ni/OneDrive/Desktop/C57_NOR/NOR_analyzed.csv'
outpath='/Users/ar_ni/OneDrive/Desktop/C57_NOR/R_Graphs/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)

df0<-subset(df, (genotype=='C57'))
df<-df0

dfLP<-data.frame(df$animal, df$genotype, df$sex, df$agegroup, df$LP)
dfLP<-cbind(dfLP, Trial='LocationPreference')
colnames(dfLP)<-c('AnimalID', 'Genotype', 'Sex', 'AgeGroup', 'Recognition.Index', 'Trial')

dfRIT2<-data.frame(df$animal, df$genotype, df$sex, df$agegroup, df$RI_T2)
dfRIT2<-cbind(dfRIT2, Trial='Trial2')
colnames(dfRIT2)<-c('AnimalID', 'Genotype', 'Sex', 'AgeGroup', 'Recognition.Index', 'Trial')

dfRIT3<-data.frame(df$animal, df$genotype, df$sex, df$agegroup, df$RI_T3)
dfRIT3<-cbind(dfRIT3, Trial='Trial3')
colnames(dfRIT3)<-c('AnimalID', 'Genotype', 'Sex', 'AgeGroup', 'Recognition.Index', 'Trial')

dfAll=rbind(dfLP, dfRIT2, dfRIT3)

#Plot mean +SSE separated by sex
ggerrorplot(dfAll, x='Trial', y='Recognition.Index', color='AgeGroup', fill='AgeGroup',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='', ylab='Recognition Index',
            legend='top', position=position_dodge(0.2),facet.by='Sex')

ggsave(paste(outpath,'NORlineC57sex.pdf', sep=''), plot=last_plot(), device='pdf', scale=1, width=5,
       height=5, units=c("in"), dpi=300)


#Plot mean +SSE
ggerrorplot(dfAll, x='Trial', y='Recognition.Index', color='AgeGroup', fill='AgeGroup',
            desc_stat='mean_se', palette=c('blueviolet', 'chartreuse1'), size=1,
            error.plot='errorbar', add='mean', point.size=1.5, xlab='', ylab='Recognition Index',
            legend='top', position=position_dodge(0.2))

ggsave(paste(outpath,'NORlineC57.pdf', sep=''), plot=last_plot(), device='pdf', scale=1, width=5,
       height=5, units=c("in"), dpi=300)

temp <- subset(dfAll, Trial %in% "Trial3")

testMethod<-t.test(Recognition.Index ~ AgeGroup, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'RIC57Trial3stats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAll, AgeGroup %in% "Older")
temp <- subset(temp, Trial %in% c("LocationPreference", "Trial2"))

testMethod<-t.test(Recognition.Index ~ Trial, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'RIC57LPTrial2Olderstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)

temp <- subset(dfAll, AgeGroup %in% "Younger")
temp <- subset(temp, Trial %in% c("Trial2", "Trial3"))

testMethod<-t.test(Recognition.Index ~ Trial, data = temp)

mytTable<-as_tibble(
  cbind(testMethod$data.name, testMethod$statistic, testMethod$p.value, testMethod$parameter[1], nrow(temp))
)

mycolnames<-c('contrast', 'statistic', 'p.value', 'df', 'observations')

myfile<-paste(outpath,'RIC57Trial2&3Youngerstats.csv')
write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
