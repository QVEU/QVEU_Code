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


ET_stat<-ET_merge%>%
  group_by(Sample) %>%
  t_test(value ~ Time) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
ET_stat


  
  stat.test <- mydata.long %>%
  group_by(Sample) %>%
  t_test(value ~ Species) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
  
ET_merge%>%group_by()%>%summarise(avg=mean("WT"))
