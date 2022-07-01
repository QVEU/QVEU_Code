library(tidyverse)
install.packages("ggbump")
library(ggbump)
setwd("~/Downloads/")

inputData<-read_csv("~/GitHub/QVEU_Code/RCode/RJW/04_12_22_HCL3_completedata.csv")

inputData<-inputData %>% group_by(construct) %>% mutate(ctrlMean1= mean(Ratio1[`[PTB^-1]`==0])) %>% mutate(normVal1=Ratio1/ctrlMean1)
inputData<-inputData %>% group_by(construct) %>% mutate(ctrlMean2= mean(Ratio2[`[PTB^-1]`==0])) %>% mutate(normVal2=Ratio2/ctrlMean2)

# check correlation of Biological Reps

#fit model using lm() - "linear models"
modelP<-summary(with(inputData,lm(Ratio2~Ratio1)))

#plot scatter plot colored by concentration of protein
ggplot(inputData,aes(Ratio1, Ratio2)) +
  geom_point(aes(col=log10(`[PTB^-1]`)))+
  geom_smooth(method = "lm")+
  ggtitle(label = paste("Biological Reps, Rsq=",round(modelP$adj.r.squared,3)))+
  scale_x_log10()+scale_y_log10()

#plot ratio as a function of concentration
ggplot(inputData)+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio2))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio2,group=construct,col=construct), method="glm")+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio1))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio1,group=construct,col=construct), method="glm",lty=2)+
  scale_x_log10()+scale_y_log10()+xlab("[PTB] (mM)")+ylab("Normalized IRES Activity\nFirefly/Renilla")

#logY Relative Activity
ggplot(inputData)+
  theme_bw()+
  geom_hline(yintercept = 100)+
  geom_text(aes(1,100,label="(No PTB)"),nudge_y = .065,size=4)+
#  geom_point(aes(`[PTB^-1]`*2.9, normVal2,col=construct))+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal2*100,col=construct),fun.y = mean,geom = "line", lwd =1)+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal2*100,col=construct),fun.data = mean_se,geom = "pointrange", lwd =1,alpha =0.4)+
  #geom_smooth(aes(`[PTB^-1]`*2.9, normVal2,group=construct,col=construct), method="glm")+
#  geom_point(aes(`[PTB^-1]`*2.9, normVal1,col=construct))+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal1*100,col=construct),fun.y = mean,geom = "line",lty=2,lwd=1)+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal1*100,col=construct),fun.data = mean_se,geom = "pointrange",lty=1,lwd=1,alpha=0.4)+
  #geom_smooth(aes(`[PTB^-1]`*2.9, normVal1,group=construct,col=construct), method="glm",lty=2)+
  scale_x_log10()+scale_y_log10()+xlab("[PTB] (mM)")+ylab("Normalized IRES Activity (%)")+scale_color_brewer(palette = "Dark2")

#linearY Relative Activity
ggplot(inputData)+
  theme_bw()+
  geom_hline(yintercept = 100)+
  geom_text(aes(1,100,label="(No PTB)"),nudge_y = 8,size=4)+
  #  geom_point(aes(`[PTB^-1]`*2.9, normVal2,col=construct))+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal2*100,col=construct),fun.y = mean,geom = "line", lwd =1)+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal2*100,col=construct),fun.data = mean_se,geom = "pointrange", lwd =1,alpha =0.4)+
  #geom_smooth(aes(`[PTB^-1]`*2.9, normVal2,group=construct,col=construct), method="glm")+
  #  geom_point(aes(`[PTB^-1]`*2.9, normVal1,col=construct))+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal1*100,col=construct),fun.y = mean,geom = "line",lty=2,lwd=1)+
  stat_summary(aes(`[PTB^-1]`*2.9,normVal1*100,col=construct),fun.data = mean_se,geom = "pointrange",lty=1,lwd=1,alpha=0.4)+
  #geom_smooth(aes(`[PTB^-1]`*2.9, normVal1,group=construct,col=construct), method="glm",lty=2)+
  scale_x_log10()+xlab("[PTB] (mM)")+ylab("Normalized IRES Activity (%)")+scale_color_brewer(type = 'div',palette = )

ggplot(inputData)+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio2/mean(inputData[construct=="negcontrol",Ratio2])))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio2,group=construct,col=construct))+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio1))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio1,group=construct,col=construct),lty=2)+
  scale_x_log10()+scale_y_log10()+xlab("[PTB] (mM)")+ylab("Normalized IRES Activity\nFirefly/Renilla")

ggplot(inputData)+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio2))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio2,group=construct,col=construct))+
  geom_point(aes(`[PTB^-1]`*2.9, Ratio1))+
  geom_smooth(aes(`[PTB^-1]`*2.9, Ratio1,group=construct,col=construct),lty=2)+
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

