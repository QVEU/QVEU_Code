library(data.table)
library(ggpubr)


luc<-fread("~/Documents/Luciferase Experiments/03_28_22_HelaCellLysate_exp1_byProt.csv")
luc=luc[,1:9]#remove analysis columns

luc[,meanFR:=mean(FF_Renilla),by=list(Sample1)]

luc

luc[,sdFR:=sd(FF_Renilla), by=list(Sample1)]

luc

#T-tests
t.test(luc[Sample1=="UORF",FF_Renilla],luc[Sample1=="del-uORF",FF_Renilla])
t.test(luc[Sample1=="pvORF",FF_Renilla],luc[Sample1=="del-pvORF",FF_Renilla])


ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Sample1,meanFR,mean),fill=Sample1),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Sample1))+
  scale_fill_discrete()+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

#Patrick
ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Genotype,meanFR,mean),fill=Genotype),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Genotype))+
  scale_fill_discrete()+
  facet_wrap(~Protein,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Protein,meanFR,mean),fill=Protein),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Protein))+
  scale_fill_discrete()+
  facet_wrap(~Genotype,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))


#ggplot(deluORF_uORF[Rep=="1A"])+
#geom_bar(aes(y=deluORF_uORF,x=Sample1),stat='identity')


#uORF and pvORF Proportion - 1
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF<-cov(luc[Sample1=="del-uORF",FF_Renilla],luc[Sample1=="UORF",FF_Renilla])
covpvORF<-cov(luc[Sample1=="del-pvORF",FF_Renilla],luc[Sample1=="pvORF",FF_Renilla])

delEffect<-data.table(comparison=c("uORF","pvORF"),
                      proportion=c(meansLuc[Sample1=="del-uORF",meanFR]/meansLuc[Sample1=="UORF",meanFR],meansLuc[Sample1=="del-pvORF",meanFR]/meansLuc[Sample1=="pvORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uORF",sdFR]/meansLuc[Sample1=="del-uORF",meanFR])**2)+((meansLuc[Sample1=="UORF",sdFR]/meansLuc[Sample1=="UORF",meanFR])**2)+2*(covUORF/(meansLuc[Sample1=="del-uORF",meanFR]*meansLuc[Sample1=="UORF",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvORF",sdFR]/meansLuc[Sample1=="del-pvORF",meanFR])**2)+((meansLuc[Sample1=="pvORF",sdFR]/meansLuc[Sample1=="pvORF",meanFR])**2)+2*(covpvORF/(meansLuc[Sample1=="del-pvORF",meanFR]*meansLuc[Sample1=="pvORF",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

##uORF-PTB1 and pvORF-PTB1 Proportion - 1
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb1<-cov(luc[Sample1=="del-uorf_ptb1",FF_Renilla],luc[Sample1=="uorf_ptb1",FF_Renilla])
covpvORF_ptb1<-cov(luc[Sample1=="del-pvorf_ptb1",FF_Renilla],luc[Sample1=="pvorf_ptb1",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_ptb1","pvorf_ptb1"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb1",meanFR]/meansLuc[Sample1=="uorf_ptb1",meanFR],meansLuc[Sample1=="del-pvorf_ptb1",meanFR]/meansLuc[Sample1=="pvorf_ptb1",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb1",sdFR]/meansLuc[Sample1=="del-uorf_ptb1",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb1",sdFR]/meansLuc[Sample1=="uorf_ptb1",meanFR])**2)+2*(covUORF_ptb1/(meansLuc[Sample1=="del-uorf_ptb1",meanFR]*meansLuc[Sample1=="uorf_ptb1",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_ptb1",sdFR]/meansLuc[Sample1=="del-pvorf_ptb1",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb1",sdFR]/meansLuc[Sample1=="pvorf_ptb1",meanFR])**2)+2*(covpvORF_ptb1/(meansLuc[Sample1=="del-pvorf_ptb1",meanFR]*meansLuc[Sample1=="pvorf_ptb1",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

##uORF-PTB2 and pvORF-PTB2 Proportion - 1
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb2<-cov(luc[Sample1=="del-uorf_ptb2",FF_Renilla],luc[Sample1=="uorf_ptb2",FF_Renilla])
covpvORF_ptb2<-cov(luc[Sample1=="del-pvorf_ptb2",FF_Renilla],luc[Sample1=="pvorf_ptb2",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_ptb2","pvorf_ptb2"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb2",meanFR]/meansLuc[Sample1=="uorf_ptb2",meanFR],meansLuc[Sample1=="del-pvorf_ptb2",meanFR]/meansLuc[Sample1=="pvorf_ptb2",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb2",sdFR]/meansLuc[Sample1=="del-uorf_ptb2",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb2",sdFR]/meansLuc[Sample1=="uorf_ptb2",meanFR])**2)+2*(covUORF_ptb2/(meansLuc[Sample1=="del-uorf_ptb2",meanFR]*meansLuc[Sample1=="uorf_ptb2",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_ptb2",sdFR]/meansLuc[Sample1=="del-pvorf_ptb2",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb2",sdFR]/meansLuc[Sample1=="pvorf_ptb2",meanFR])**2)+2*(covpvORF_ptb2/(meansLuc[Sample1=="del-pvorf_ptb2",meanFR]*meansLuc[Sample1=="pvorf_ptb2",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


##uORF-PTB4 and pvORF-PTB4 Proportion - 1
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb4<-cov(luc[Sample1=="del-uorf_ptb4",FF_Renilla],luc[Sample1=="uorf_ptb4",FF_Renilla])
covpvORF_ptb4<-cov(luc[Sample1=="del-pvorf_ptb4",FF_Renilla],luc[Sample1=="pvorf_ptb4",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_ptb4","pvorf_ptb4"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb4",meanFR]/meansLuc[Sample1=="uorf_ptb4",meanFR],meansLuc[Sample1=="del-pvorf_ptb4",meanFR]/meansLuc[Sample1=="pvorf_ptb4",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb4",sdFR]/meansLuc[Sample1=="del-uorf_ptb4",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb4",sdFR]/meansLuc[Sample1=="uorf_ptb4",meanFR])**2)+2*(covUORF_ptb4/(meansLuc[Sample1=="del-uorf_ptb4",meanFR]*meansLuc[Sample1=="uorf_ptb4",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_ptb4",sdFR]/meansLuc[Sample1=="del-pvorf_ptb4",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb4",sdFR]/meansLuc[Sample1=="pvorf_ptb4",meanFR])**2)+2*(covpvORF_ptb4/(meansLuc[Sample1=="del-pvorf_ptb4",meanFR]*meansLuc[Sample1=="pvorf_ptb4",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#uORF-PTB4 and pvORF-PTB4 Proportion
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb4<-cov(luc[Sample1=="del-uorf_ptb4",FF_Renilla],luc[Sample1=="uorf_ptb4",FF_Renilla])
covpvORF_ptb4<-cov(luc[Sample1=="del-pvorf_ptb4",FF_Renilla],luc[Sample1=="pvorf_ptb4",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_ptb4","pvorf_ptb4"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb4",meanFR]/meansLuc[Sample1=="uorf_ptb4",meanFR],meansLuc[Sample1=="del-pvorf_ptb4",meanFR]/meansLuc[Sample1=="pvorf_ptb4",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb4",sdFR]/meansLuc[Sample1=="del-uorf_ptb4",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb4",sdFR]/meansLuc[Sample1=="uorf_ptb4",meanFR])**2)+2*(covUORF_ptb4/(meansLuc[Sample1=="del-uorf_ptb4",meanFR]*meansLuc[Sample1=="uorf_ptb4",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_ptb4",sdFR]/meansLuc[Sample1=="del-pvorf_ptb4",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb4",sdFR]/meansLuc[Sample1=="pvorf_ptb4",meanFR])**2)+2*(covpvORF_ptb4/(meansLuc[Sample1=="del-pvorf_ptb4",meanFR]*meansLuc[Sample1=="pvorf_ptb4",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()
