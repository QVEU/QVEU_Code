library(tidyverse)
library(dplyr)
library(rstatix)
library(ggpubr)

ET_Time<-read_delim("/Volumes/lvd-qveu/Rebecca Walker/7-1-22-NEO2_Replicon_timecourse/Renilla_EnduRen_Timecourse_edited.txt")
ET_long<-tidyr::pivot_longer(ET_Time,cols = 3:98)
ET_long[,c("Row","Col")]<-str_split_fixed(ET_long$name, "", 2)

labels<-tibble(Col=c(1:12),
               Sample=c("water","WT","WT","del564-568","del564-568","WT","WT","del564-568","del564-568","mock","mock","water"),
               Trx=c("water","none","none","none","none","guan","guan","guan","guan","none","none","water"))

ET_merge<-merge(ET_long,labels,by = "Col")
ET_merge$G=factor(paste(ET_merge$Sample, ET_merge$Trx))
ET_merge = ET_merge%>%group_by(Time,G) %>%mutate(means=mean(value),sem=sd(value)/sqrt(length(value)),len=(length(value)))

ET_stat<-ET_merge%>%ungroup()%>%filter(G%in%c("WT none","del564-568 none"))%>%
  group_by(Time) %>%
  t_test(value ~ G) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance()
ET_stat

ET_stat_guan<-ET_merge%>%ungroup()%>%filter(G%in%c("WT guan","del564-568 guan"))%>%
  group_by(Time) %>%
  t_test(value ~ G) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance()



  
ET_merge%>%group_by()%>%summarise(avg=mean(Sample))
