neo2data<-read.csv("/Volumes/lvd-qveu/Rebecca_Walker/neo2/22_12_20/22_12_20_e144a_rep3.csv")
samplesheet<-read.csv("/Volumes/lvd-qveu/Rebecca_Walker/neo2/22_12_20/22_12_20_samplesheet.csv")

library(tidyr)
library(datetime)
library(stringr)
library(dplyr)
library(ggplot2)

neo2data

neo2_long<-pivot_longer(data=neo2data,cols=seq(3,98),names_to="well", values_to="luminescence")
neo2_merge<-merge(samplesheet,neo2_long, by="well")
neo2_merge

ggplot(neo2_merge)+geom_point(aes(as.time(Time),luminescence))

ggplot(neo2_merge)+
  #geom_point(aes(as.time(Time),luminescence,col=sample))+
  geom_line(aes(as.time(Time),luminescence,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time),luminescence, group=sample),lwd=2)+
  scale_color_discrete()#+scale_y_log10()

ggplot(neo2_merge)+
  #geom_point(aes(as.time(Time),luminescence,col=sample))+
  geom_line(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
  xlab("Time (hours)")+
  scale_color_discrete()#+scale_y_log10()

ggplot(
  neo2_merge[neo2_merge$sample%in%c("del E144A Renilla + WT Renilla + Guan","del Renilla + Guan","WT Renilla + Guan","WT E144A Renilla + del Renilla + Guan","mock","del E144A Renilla","WT E144A Renilla"),]
  )+
  #geom_point(aes(as.time(Time),luminescence,col=sample))+
  geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
  xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()


ggplot(
    neo2_merge[neo2_merge$Time%in%c("0:00:00",	
                                    "0:30:00",
                                    "1:00:00",	
                                    "1:30:00",	
                                    "2:00:00",	
                                    "2:30:00",	
                                    "3:00:00",	
                                    "3:30:00",	
                                    "4:00:00",	
                                    "4:30:00",	
                                    "5:00:00",	
                                    "5:30:00",
                                    "6:00:00",	
                                    "6:30:00",	
                                    "7:00:00",	
                                    "7:30:00",	
                                    "8:00:00"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  
  
  
  
ggplot(
    neo2_merge[neo2_merge$Time%in%c("0:00:00",	
                                    "0:30:00",
                                    "1:00:00",	
                                    "1:30:00",	
                                    "2:00:00",	
                                    "2:30:00",	
                                    "3:00:00",	
                                    "3:30:00",	
                                    "4:00:00",	
                                    "4:30:00",	
                                    "5:00:00",	
                                    "5:30:00",
                                    "6:00:00",	
                                    "6:30:00"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  
  
  
ggplot(
    neo2_merge[neo2_merge$Time%in%c("0:00:00",	
                                    "0:30:00",
                                    "1:00:00",	
                                    "1:30:00",	
                                    "2:00:00",	
                                    "2:30:00",	
                                    "3:00:00",	
                                    "3:30:00",	
                                    "4:00:00",	
                                    "4:30:00",	
                                    "5:00:00"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()  
  
  
ggplot(
    neo2_merge[neo2_merge$sample%in%c("del E144A Renilla + WT Renilla","WT E144A Renilla + del Renilla"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  
ggplot(
    neo2_merge[neo2_merge$sample%in%c("WT Renilla","del Renilla"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  
ggplot(
    neo2_merge[neo2_merge$sample%in%c("del E144A Renilla + WT Renilla","WT E144A Renilla + del Renilla","WT Renilla","del Renilla","WT Renilla (2x RNA)"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  
  ggplot(
    neo2_merge[(as.time(neo2_merge$Time)<32500)&neo2_merge$sample%in%c("del E144A Renilla + WT Renilla","WT E144A Renilla + del Renilla","WT Renilla","del Renilla","WT Renilla (2x RNA)"),]
  )+
    #geom_point(aes(as.time(Time),luminescence,col=sample))+
    geom_point(aes(as.time(Time)/3600,luminescence,col=sample,group=well),alpha=0.4)+
    stat_summary(geom = "line",fun = mean,aes(col=sample,as.time(Time)/3600,luminescence, group=sample),lwd=2)+
    xlab("Time (hours)")
  scale_color_discrete()#+scale_y_log10()
  neo2_merge$inttime<-as.integer(as.time(neo2_merge$Time))  
  write.csv(neo2_merge,quote = F,row.names = F,file = "annotatedData_del564-8_replicon_neo2.csv")

  
neo2_TimeCorrected<-neo2_merge(as.time(Time)/3600)



