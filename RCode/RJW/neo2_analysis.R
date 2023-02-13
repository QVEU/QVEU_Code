
library(tidyr)
library(datetime)
library(stringr)
library(dplyr)
library(ggplot2)

#read in csv file and convert time data to seconds
neo2_merge <-read.csv("/Volumes/lvd-qveu/Rebecca_Walker/neo2/Merged_Neo2_Files.csv")
neo2_merge$Time_sec<-as.integer(as.time(neo2_merge$Time))


ggplot(neo2_merge[neo2_merge$sample!="negative",])+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge)+
  facet_grid(cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge)+
  facet_grid(cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+xlab("Time (hours)")

ggplot(neo2_merge)+
  facet_grid(rows = vars(sample),cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sample%in%c("del E144A Renilla + WT Renilla","WT E144A Renilla + del Renilla"),])+
  facet_grid(rows = vars(sample),cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sample%in%c("WT E144A Renilla + del Renilla GDD mutpol","WT E144A Renilla + WT Renilla GDD mutpol",
                                         "del E144A Renilla + del Renilla GDD mutpol","del E144A Renilla + WT Renilla GDD mutpol"),])+
  facet_grid(rows = vars(sample),cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sample%in%c("WT E144A Renilla + del Renilla GDD mutpol","WT E144A Renilla + WT Renilla GDD mutpol",
                                         "del E144A Renilla + del Renilla GDD mutpol","del E144A Renilla + WT Renilla GDD mutpol"),])+
  facet_grid(rows = vars(sample),cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sample%in%c("WT E144A Renilla + del Renilla GDD mutpol","WT E144A Renilla + WT Renilla GDD mutpol",
                                         "del E144A Renilla + del Renilla GDD mutpol","del E144A Renilla + WT Renilla GDD mutpol"),])+
  facet_grid(cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+xlim(0,9)+scale_y_log10()+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sample%in%c("WT E144A Renilla + del Renilla GDD mutpol","WT E144A Renilla + WT Renilla GDD mutpol",
                                         "del E144A Renilla + del Renilla GDD mutpol","del E144A Renilla + WT Renilla GDD mutpol"),])+
  facet_grid(cols = vars(sampleID))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+xlim(0,9)+xlab("Time (hours)")

ggplot(neo2_merge[neo2_merge$sampleID%in%c("23_02_09_renillaluc_replicon.txt","23_02_10_renilla_replicon.txt"),])+
  facet_grid(rows=vars(sampleID),cols = vars(sample))+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count, group=well),lwd=2)+
  scale_color_discrete()+scale_y_log10()+xlab("Time (hours)")


 # I am confused
ggplot(neo2_merge[neo2_merge$sampleID%in%c("23_02_09_renillaluc_replicon.txt"),][neo2_merge$sample%in%c("del Renilla","WT Renilla"),])+
  geom_line(aes((Time_sec-TimeDelay_sec)/3600,count,col=sample,group=well),alpha=0.4)+
  stat_summary(geom = "line",fun = mean,aes(col=sample,(Time_sec-TimeDelay_sec)/3600,count,group=well),lwd=2)+
  scale_color_discrete()+xlim(0,8)+scale_y_log10()+xlab("Time (hours)")
  
