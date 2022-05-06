library(data.table)
library(reshape2)
library(limma)
library(ggpubr)

info<-fread("/Users/walkerrj/Documents/Luciferase Experiments/04_12_22_HCL3_Data/04_12_22_HCL3_table.txt",header=TRUE,stringsAsFactors = TRUE)
info$variable<-as.factor(info$variable)
info$RepB<-as.factor(info$RepB)
alldata=data.table()
for(file in list.files(full.names = TRUE, pattern=".csv","/Users/walkerrj/Documents/Luciferase\ Experiments/04_12_22_HCL3_Data")){
  f=fread(header = TRUE,file,stringsAsFactors = TRUE)
  print(f)
  m<-(melt(f))
  m$datatable=file
  alldata=rbindlist(list(alldata,m))
  }
alldata[,reading:=strsplit2(strsplit2(datatable,"/Users/walkerrj/Documents/Luciferase Experiments/04_12_22_HCL3_Data/04_12_22_HCL3_")[,2],".csv")[,1]]
alldata[,RepB:=strsplit2(reading,"_")[,2]]

mergedInfo<-merge.data.table(alldata,info,by = c("V1","variable","RepB"))

mergedInfo<-mergedInfo[,reading:=strsplit2(strsplit2(datatable,"/Users/walkerrj/Documents/Luciferase Experiments/04_12_22_HCL3_Data/04_12_22_HCL3_")[,2],".csv")[,1]]

castdata<-dcast.data.table(mergedInfo[,-c("RepB","datatable")],value.var = "value",formula = ...~reading)
DT<-data.table(castdata[,7:10]-colMeans(castdata[construct=="negcontrol",7:10]))
colnames(DT)<-c("FF1_adj","FF2_adj","R1_adj","R2_adj")
castdata<-cbind(castdata,DT)

castdata[,Ratio1:=FF1_adj/R1_adj]
castdata[,Ratio2:=FF2_adj/R2_adj]
castdata[,sdfr1:=sd(Ratio1),by=c("construct","`[PTB^-1]`")]
castdata[,sdfr2:=sd(Ratio2),by=c("construct","`[PTB^-1]`")]
castdata[,meanfr1:=mean(Ratio1),by=c("construct","`[PTB^-1]`")]
castdata[,meanfr2:=mean(Ratio2),by=c("construct","`[PTB^-1]`")]


castdata

library(ggplot2)
ggplot(castdata)+geom_line(aes(`[PTB^-1]`,Ratio1,group=RepT),col="cyan")+
  geom_line(aes(`[PTB^-1]`,Ratio2,group=RepT),col="red")+
  facet_wrap(~construct,scales = "free_y")+scale_x_log10()+ylab("FF:Renilla")

ggplot(castdata)+geom_line(aes(`[PTB^-1]`,FF1_adj,group=RepT),col="cyan")+
  geom_line(aes(`[PTB^-1]`,FF2_adj,group=RepT),col="red")+
  facet_wrap(~construct,scales = "free_y")+scale_x_log10()+ylab("FF")

ggplot(castdata)+geom_line(aes(`[PTB^-1]`,R1_adj,group=RepT),col="cyan")+
  geom_line(aes(`[PTB^-1]`,R2_adj,group=RepT),col="red")+
  facet_wrap(~construct,scales = "free_y")+scale_x_log10()+ylab("FF")

#working large deletion effect graph

meansdata<-castdata[,list(meanfr1,sdfr1),by=c("construct","`[PTB^-1]`")]
covUORF<-cov(castdata[construct=="del-uorf",Ratio1],castdata[construct=="uorf",Ratio1])
covpvORF<-cov(castdata[construct=="del-pvorf",Ratio1],castdata[construct=="pvorf",Ratio1])

delEffect<-data.table(comparison=c("uorf","pvorf"),dilution=c("[PTB^-1]"),
                      proportion=c(meansdata[construct=="del-uorf",meanfr1]/meansdata[construct=="uorf",meanfr1],meansdata[construct=="del-pvorf",meanfr1]/meansdata[construct=="pvorf",meanfr1]),
                      Error=c(sqrt((meansdata[construct=="del-uorf",sdfr1]/meansdata[construct=="del-uorf",meanfr1])**2)+((meansdata[construct=="uorf",sdfr1]/meansdata[construct=="uorf",meanfr1])**2)+2*(covUORF/(meansdata[construct=="del-uorf",meanfr1]*meansdata[construct=="uorf",meanfr1])),
                              sqrt((meansdata[construct=="del-pvorf",sdfr1]/meansdata[construct=="del-pvorf",meanfr1])**2)+((meansdata[construct=="pvorf",sdfr1]/meansdata[construct=="pvorf",meanfr1])**2)+2*(covpvORF/(meansdata[construct=="del-pvorf",meanfr1]*meansdata[construct=="pvorf",meanfr1]))
                      ))

ggplot(delEffect)+geom_line(aes(y=proportion-1,x=comparison))

ggplot(delEffect, aes(y=proportion-1,x=comparison, colour = "construct"),stat='identity') + geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  geom_point()+ scale_x_discrete()

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()



#working deletion effect comparison

#deletion effect pvorf and uorf at [1] ptb
meansLuc<-castdata[`[PTB]`=="1",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="1",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="1",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="1",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="1",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#deletion effect pvorf and uorf at [4] ptb
meansLuc<-castdata[`[PTB]`=="4",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="4",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="4",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="4",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="4",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#deletion effect pvorf and uorf at [16] ptb
meansLuc<-castdata[`[PTB]`=="16",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="16",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="16",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="16",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="16",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#deletion effect pvorf and uorf at [64] ptb
meansLuc<-castdata[`[PTB]`=="64",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="64",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="64",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="64",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="64",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#deletion effect pvorf and uorf at [256] ptb
meansLuc<-castdata[`[PTB]`=="256",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="256",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="256",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="256",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="256",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#deletion effect pvorf and uorf at [1024] ptb
meansLuc<-castdata[`[PTB]`=="1024",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="1024",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="1024",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="1024",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="1024",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#deletion effect pvorf and uorf at [0] ptb
meansLuc<-castdata[`[PTB]`=="0",list(meanFR=mean(Ratio1),sdFR=sd(Ratio1)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="0",Ratio1],castdata[construct=="del-pvorf"&`[PTB]`=="0",Ratio1])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="0",Ratio1],castdata[construct=="del-uorf"&`[PTB]`=="0",Ratio1])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#deletion effect pvorf and uorf at [0] ptb
meansLuc<-castdata[`[PTB]`=="0",list(meanFR=mean(Ratio1&Ratio2),sdFR=sd(Ratio1&Ratio2)),by=list(construct)]
covpvORF<-cov(castdata[construct=="pvorf"&`[PTB]`=="0",Ratio1,Ratio2],castdata[construct=="del-pvorf"&`[PTB]`=="0",Ratio1,Ratio2])
covuORF<-cov(castdata[construct=="uorf"&`[PTB]`=="0",Ratio1,Ratio2],castdata[construct=="del-uorf"&`[PTB]`=="0",Ratio1,Ratio2])

delEffect<-data.table(comparison=c("pvorf","uorf"),
                      proportion=c(meansLuc[construct=="del-uorf",meanFR]/meansLuc[construct=="uorf",meanFR],meansLuc[construct=="del-pvorf",meanFR]/meansLuc[construct=="pvorf",meanFR]),
                      Error=c(sqrt((meansLuc[construct=="del-uorf",sdFR]/meansLuc[construct=="del-uorf",meanFR])**2)+((meansLuc[construct=="uorf",sdFR]/meansLuc[construct=="uorf",meanFR])**2)+2*(covuORF/(meansLuc[construct=="del-pvorf",meanFR]*meansLuc[construct=="uorf",meanFR])),
                              sqrt((meansLuc[construct=="del-pvorf",sdFR]/meansLuc[construct=="del-pvorf",meanFR])**2)+((meansLuc[construct=="pvorf",sdFR]/meansLuc[construct=="pvorf",meanFR])**2)+2*(covpvORF/(meansLuc[construct=="del-pvORF",meanFR]*meansLuc[construct=="pvorf",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#Working on figure RW
ggplot(castdata, aes(`[PTB^-1]`,Ratio1, colour = construct)) +
  geom_point()

ggplot(castdata, aes(`[PTB^-1]`,Ratio2, colour = construct)) +
  geom_point()
