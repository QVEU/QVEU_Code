library(tidyverse)
#install.packages("ggbump")
library(ggbump)
setwd("~/Downloads/")
inputData<-read_csv("./04_12_22_HCL3_completedata.csv")

# check correlation of Biological Reps

#fit model using lm() - "linear models"
modelP<-summary(with(inputData,lm(Ratio2~Ratio1)))

#plot scatter plot colored by concentration of protein
ggplot(inputData,aes(Ratio1, Ratio2))+
  geom_point(aes(col=log10(`[PTB^-1]`)))+
  geom_smooth(method = "lm")+
  ggtitle(label = paste("Biological Reps, Rsq=",round(modelP$adj.r.squared,3)))+
  scale_x_log10()+scale_y_log10()

#plot ratio as a function of concentration
ggplot(inputData)+
  geom_point(aes(`[PTB^-1]`*5, Ratio2))+
  geom_smooth(aes(`[PTB^-1]`*5, Ratio2,group=construct,col=construct), method="glm")+
  geom_point(aes(`[PTB^-1]`*5, Ratio1))+
  geom_smooth(aes(`[PTB^-1]`*5, Ratio1,group=construct,col=construct), method="glm",lty=2)+
  scale_x_log10()+scale_y_log10()+xlab("[PTB] (mM)")+ylab("Normalized IRES Activity\nFirefly/Renilla")

ggplot(inputData)+
  geom_bar(aes(`[PTB^-1]`+0.0001, FF1_adj),stat="identity", fill="lightgreen",position='dodge')+scale_x_log10()+
  geom_bar(aes(`[PTB^-1]`+0.0001, R1_adj),stat="identity", fill="darkgreen",position='dodge')+facet_wrap(~construct,scales='free')
  geom_smooth(aes(`[PTB^-1]`, FF1_adj,group=construct,col=construct), method="glm")+
  geom_point(aes(`[PTB^-1]`, FF2_adj))+
  geom_smooth(aes(`[PTB^-1]`, FF2_adj ,group=construct,col=construct), method="glm")+
  scale_y_log10()+scale_x_log10()

ggplot(inputData)+
  geom_point(aes(`[PTB^-1]`, R1_adj))+
  geom_smooth(aes(`[PTB^-1]`, R1_adj,group=construct,col=construct), method="glm")+
  geom_point(aes(`[PTB^-1]`, R2_adj))+
  geom_smooth(aes(`[PTB^-1]`, R2_adj ,group=construct,col=construct), method="glm")+
  scale_y_log10()+scale_x_log10()

ggplot(inputData,aes(Ratio1, Ratio2))+
  geom_point()+
  geom_smooth(method = "lm")+
  ggtitle(label = paste("Biological Reps, Rsq=",modelP$adj.r.squared))+
  scale_x_log10()+scale_y_log10()

