library(data.table)
library(ggpubr)


luc<-fread("~/Documents/Luciferase Experiments/03_31_22_HelaCellLysate_exp2_byProt.csv")
luc=luc[,1:10]#remove analysis columns

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

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Accesory,meanFR,mean),fill=Accesory),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Accesory))+
  scale_fill_discrete()+
  facet_wrap(~Genotype,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Genotype,meanFR,mean),fill=Accesory),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Genotype))+
  scale_fill_discrete()+
  facet_wrap(~Protein,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Sample1,meanFR,mean),fill=Genotype),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Sample1))+
  scale_fill_discrete()+
  facet_wrap(~Accesory,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Sample1,meanFR,mean),fill=Genotype),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Sample1))+
  scale_fill_discrete()+
  facet_wrap(~Protein,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Sample1,meanFR,mean),fill=Protein),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Sample1))+
  scale_fill_discrete()+
  facet_wrap(~Genotype,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

ggplot(luc[Sample1!="control-"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Genotype,meanFR,mean),fill=Protein),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Genotype))+
  scale_fill_discrete()+
  facet_wrap(~Protein,scales="free")+
  #scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")+theme(axis.text.x = element_text(hjust = 1, angle=50))

#ggplot(deluORF_uORF[Rep=="1A"])+
#geom_bar(aes(y=deluORF_uORF,x=Sample1),stat='identity')

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

#PTB2
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

#PTB2 and accessory proteins
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb2_acc<-cov(luc[Sample1=="del-uorf_ptb2_acc",FF_Renilla],luc[Sample1=="uorf_ptb2_acc",FF_Renilla])
covpvORF_ptb2_acc<-cov(luc[Sample1=="del-pvorf_ptb2_acc",FF_Renilla],luc[Sample1=="pvorf_ptb2_acc",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_ptb2_acc","pvorf_ptb2_acc"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb2_acc",meanFR]/meansLuc[Sample1=="uorf_ptb2_acc",meanFR],meansLuc[Sample1=="del-pvorf_ptb2_acc",meanFR]/meansLuc[Sample1=="pvorf_ptb2_acc",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb2_acc",sdFR]/meansLuc[Sample1=="del-uorf_ptb2_acc",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb2_acc",sdFR]/meansLuc[Sample1=="uorf_ptb2_acc",meanFR])**2)+2*(covUORF_ptb2_acc/(meansLuc[Sample1=="del-uorf_ptb2_acc",meanFR]*meansLuc[Sample1=="uorf_ptb2_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_ptb2_acc",sdFR]/meansLuc[Sample1=="del-pvorf_ptb2_acc",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb2_acc",sdFR]/meansLuc[Sample1=="pvorf_ptb2_acc",meanFR])**2)+2*(covpvORF_ptb2_acc/(meansLuc[Sample1=="del-pvorf_ptb2_acc",meanFR]*meansLuc[Sample1=="pvorf_ptb2_acc",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#Accesory Proteins
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_acc<-cov(luc[Sample1=="del-uorf_acc",FF_Renilla],luc[Sample1=="uorf_acc",FF_Renilla])
covpvORF_acc<-cov(luc[Sample1=="del-pvorf_acc",FF_Renilla],luc[Sample1=="pvorf_acc",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_acc","pvorf_acc"),
                      proportion=c(meansLuc[Sample1=="del-uorf_acc",meanFR]/meansLuc[Sample1=="uorf_acc",meanFR],meansLuc[Sample1=="del-pvorf_acc",meanFR]/meansLuc[Sample1=="pvorf_acc",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_acc",sdFR]/meansLuc[Sample1=="del-uorf_acc",meanFR])**2)+((meansLuc[Sample1=="uorf_acc",sdFR]/meansLuc[Sample1=="uorf_acc",meanFR])**2)+2*(covUORF_acc/(meansLuc[Sample1=="del-uorf_acc",meanFR]*meansLuc[Sample1=="uorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_acc",sdFR]/meansLuc[Sample1=="del-pvorf_acc",meanFR])**2)+((meansLuc[Sample1=="pvorf_acc",sdFR]/meansLuc[Sample1=="pvorf_acc",meanFR])**2)+2*(covpvORF_acc/(meansLuc[Sample1=="del-pvorf_acc",meanFR]*meansLuc[Sample1=="pvorf_acc",meanFR]))
                      ))

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#Additional Comparisons
#meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
#covUORF_acc_compare<-cov(luc[Sample1=="del-uorf_acc",FF_Renilla],luc[Sample1=="uorf_acc",FF_Renilla])
#covpvORF_acc_compare<-cov(luc[Sample1=="del-pvorf_acc",FF_Renilla],luc[Sample1=="pvorf_acc",FF_Renilla])
#covUORF<-cov(luc[Sample1=="del-uORF",FF_Renilla],luc[Sample1=="UORF",FF_Renilla])
#covpvORF<-cov(luc[Sample1=="del-pvORF",FF_Renilla],luc[Sample1=="pvORF",FF_Renilla])


#delEffect<-data.table(comparison=c("uorf_acc","pvorf_acc", "UORF","pvORF"),
                      proportion=c(meansLuc[Sample1=="del-uorf_acc",meanFR]/meansLuc[Sample1=="del-uORF",meanFR],meansLuc[Sample1=="del-pvorf_acc",meanFR]/meansLuc[Sample1=="del-pvORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_acc",sdFR]/meansLuc[Sample1=="del-uORF",meanFR])**2)+((meansLuc[Sample1=="uorf_acc",sdFR]/meansLuc[Sample1=="UORF",meanFR])**2)+2*(covUORF_acc/(meansLuc[Sample1=="del-uorf_acc",meanFR]*meansLuc[Sample1=="uorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvorf_acc",sdFR]/meansLuc[Sample1=="del-pvORF",meanFR])**2)+((meansLuc[Sample1=="pvorf_acc",sdFR]/meansLuc[Sample1=="pvORF",meanFR])**2)+2*(covpvORF_acc/(meansLuc[Sample1=="del-pvorf_acc",meanFR]*meansLuc[Sample1=="pvorf_acc",meanFR]))
                      ))

#ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  #scale_fill_manual(values = c("cyan", "blue","pink","orange"))+
  #theme_pubr()


#uorf vs uorf+accesory proteins
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_acc<-cov(luc[Sample1=="del-uorf_acc",FF_Renilla],luc[Sample1=="uorf_acc",FF_Renilla])
covUORF<-cov(luc[Sample1=="del-uORF",FF_Renilla],luc[Sample1=="UORF",FF_Renilla])


delEffect<-data.table(comparison=c("uorf_acc","UORF"),
                      proportion=c(meansLuc[Sample1=="del-uorf_acc",meanFR]/meansLuc[Sample1=="uorf_acc",meanFR],meansLuc[Sample1=="del-uORF",meanFR]/meansLuc[Sample1=="UORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_acc",sdFR]/meansLuc[Sample1=="del-uorf_acc",meanFR])**2)+((meansLuc[Sample1=="uorf_acc",sdFR]/meansLuc[Sample1=="uorf_acc",meanFR])**2)+2*(covUORF_acc/(meansLuc[Sample1=="del-uorf_acc",meanFR]*meansLuc[Sample1=="uorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-uORF",sdFR]/meansLuc[Sample1=="del-uORF",meanFR])**2)+((meansLuc[Sample1=="UORF",sdFR]/meansLuc[Sample1=="UORF",meanFR])**2)+2*(covUORF/(meansLuc[Sample1=="del-uORF",meanFR]*meansLuc[Sample1=="UORF",meanFR]))
                      ))
  

ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#UORF vs uorf + ptb2
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covUORF_ptb2<-cov(luc[Sample1=="del-uorf_ptb2",FF_Renilla],luc[Sample1=="uorf_ptb2",FF_Renilla])
covUORF<-cov(luc[Sample1=="del-uORF",FF_Renilla],luc[Sample1=="UORF",FF_Renilla])


delEffect<-data.table(comparison=c("uorf_ptb2","UORF"),
                      proportion=c(meansLuc[Sample1=="del-uorf_ptb2",meanFR]/meansLuc[Sample1=="uorf_ptb2",meanFR],meansLuc[Sample1=="del-uORF",meanFR]/meansLuc[Sample1=="UORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_ptb2",sdFR]/meansLuc[Sample1=="del-uorf_ptb2",meanFR])**2)+((meansLuc[Sample1=="uorf_ptb2",sdFR]/meansLuc[Sample1=="uorf_ptb2",meanFR])**2)+2*(covUORF_ptb2/(meansLuc[Sample1=="del-uorf_ptb2",meanFR]*meansLuc[Sample1=="uorf_ptb2",meanFR])),
                              sqrt((meansLuc[Sample1=="del-uORF",sdFR]/meansLuc[Sample1=="del-uORF",meanFR])**2)+((meansLuc[Sample1=="UORF",sdFR]/meansLuc[Sample1=="UORF",meanFR])**2)+2*(covUORF/(meansLuc[Sample1=="del-uORF",meanFR]*meansLuc[Sample1=="UORF",meanFR]))
                      ))


ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#pvORF vs pvorf +acc
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covpvORF_acc<-cov(luc[Sample1=="del-pvorf_acc",FF_Renilla],luc[Sample1=="pvorf_acc",FF_Renilla])
covpvORF<-cov(luc[Sample1=="del-pvORF",FF_Renilla],luc[Sample1=="pvORF",FF_Renilla])


delEffect<-data.table(comparison=c("pvorf_acc","pvORF"),
                      proportion=c(meansLuc[Sample1=="del-pvorf_acc",meanFR]/meansLuc[Sample1=="pvorf_acc",meanFR],meansLuc[Sample1=="del-pvorf",meanFR]/meansLuc[Sample1=="pvORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-pvorf_acc",sdFR]/meansLuc[Sample1=="del-pvorf_acc",meanFR])**2)+((meansLuc[Sample1=="pvorf_acc",sdFR]/meansLuc[Sample1=="pvorf_acc",meanFR])**2)+2*(covpvORF_acc/(meansLuc[Sample1=="del-pvorf_acc",meanFR]*meansLuc[Sample1=="pvorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvORF",sdFR]/meansLuc[Sample1=="del-pvORF",meanFR])**2)+((meansLuc[Sample1=="pvORF",sdFR]/meansLuc[Sample1=="pvORF",meanFR])**2)+2*(covpvORF/(meansLuc[Sample1=="del-pvORF",meanFR]*meansLuc[Sample1=="pvORF",meanFR]))
                      ))


ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#pvorf vs pvorf+ptb2
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covpvORF_ptb2<-cov(luc[Sample1=="del-pvorf_ptb2",FF_Renilla],luc[Sample1=="pvorf_ptb2",FF_Renilla])
covpvORF<-cov(luc[Sample1=="del-pvORF",FF_Renilla],luc[Sample1=="pvORF",FF_Renilla])


delEffect<-data.table(comparison=c("pvorf_ptb2","pvORF"),
                      proportion=c(meansLuc[Sample1=="del-pvorf_ptb2",meanFR]/meansLuc[Sample1=="pvorf_ptb2",meanFR],meansLuc[Sample1=="del-pvorf",meanFR]/meansLuc[Sample1=="pvORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-pvorf_ptb2",sdFR]/meansLuc[Sample1=="del-pvorf_ptb2",meanFR])**2)+((meansLuc[Sample1=="pvorf_ptb2",sdFR]/meansLuc[Sample1=="pvorf_ptb2",meanFR])**2)+2*(covpvORF_ptb2/(meansLuc[Sample1=="del-pvorf_ptb2",meanFR]*meansLuc[Sample1=="pvorf_ptb2",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvORF",sdFR]/meansLuc[Sample1=="del-pvORF",meanFR])**2)+((meansLuc[Sample1=="pvORF",sdFR]/meansLuc[Sample1=="pvORF",meanFR])**2)+2*(covpvORF/(meansLuc[Sample1=="del-pvORF",meanFR]*meansLuc[Sample1=="pvORF",meanFR]))
                      ))


ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()

#pvorf with accesorry proteins vs pvorf with accessory and ptb2
meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covpvORF_acc<-cov(luc[Sample1=="del-pvorf_acc",FF_Renilla],luc[Sample1=="pvorf_acc",FF_Renilla])
covpvORF_ptb2_acc<-cov(luc[Sample1=="del-pvorf_ptb2_acc",FF_Renilla],luc[Sample1=="pvorf_ptb2_acc",FF_Renilla])

delEffect<-data.table(comparison=c("pvorf_acc","pvORF_ptb2_acc"),
                      proportion=c(meansLuc[Sample1=="del-pvorf_acc",meanFR]/meansLuc[Sample1=="pvorf_acc",meanFR],meansLuc[Sample1=="del-pvorf_ptb2_acc",meanFR]/meansLuc[Sample1=="pvORF_ptb2_acc",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-pvorf_acc",sdFR]/meansLuc[Sample1=="del-pvorf_acc",meanFR])**2)+((meansLuc[Sample1=="pvorf_acc",sdFR]/meansLuc[Sample1=="pvorf_acc",meanFR])**2)+2*(covpvORF_acc/(meansLuc[Sample1=="del-pvorf_acc",meanFR]*meansLuc[Sample1=="pvorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-pvORF_ptb2_acc",sdFR]/meansLuc[Sample1=="del-pvORF_ptb2_acc",meanFR])**2)+((meansLuc[Sample1=="pvORF_ptb2_acc",sdFR]/meansLuc[Sample1=="pvORF_ptb2_acc",meanFR])**2)+2*(covpvORF_ptb2_acc/(meansLuc[Sample1=="del-pvORF_ptb2_acc",meanFR]*meansLuc[Sample1=="pvORF_ptb2_acc",meanFR]))
                      ))


ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#uorf with accesory vs uorf with accessory and ptb2

meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
covuorf_acc<-cov(luc[Sample1=="del-uorf_acc",FF_Renilla],luc[Sample1=="uorf_acc",FF_Renilla])
covuORF_ptb2_acc<-cov(luc[Sample1=="del-uorf_ptb2_acc",FF_Renilla],luc[Sample1=="uorf_ptb2_acc",FF_Renilla])

delEffect<-data.table(comparison=c("uorf_acc","uORF_ptb2_acc"),
                      proportion=c(meansLuc[Sample1=="del-uorf_acc",meanFR]/meansLuc[Sample1=="uorf_acc",meanFR],meansLuc[Sample1=="del-uorf_ptb2_acc",meanFR]/meansLuc[Sample1=="uORF_ptb2_acc",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_acc",sdFR]/meansLuc[Sample1=="del-uorf_acc",meanFR])**2)+((meansLuc[Sample1=="uorf_acc",sdFR]/meansLuc[Sample1=="uorf_acc",meanFR])**2)+2*(covuorf_acc/(meansLuc[Sample1=="del-uorf_acc",meanFR]*meansLuc[Sample1=="uorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-uORF_ptb2_acc",sdFR]/meansLuc[Sample1=="del-uORF_ptb2_acc",meanFR])**2)+((meansLuc[Sample1=="uORF_ptb2_acc",sdFR]/meansLuc[Sample1=="uORF_ptb2_acc",meanFR])**2)+2*(covuORF_ptb2_acc/(meansLuc[Sample1=="del-uORF_ptb2_acc",meanFR]*meansLuc[Sample1=="uORF_ptb2_acc",meanFR]))
                      ))


ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue"))+
  theme_pubr()


#comparison fix
#meansLuc<-luc[,list(meanFR=mean(FF_Renilla),sdFR=sd(FF_Renilla)),by=Sample1]
#covUORF_acc<-cov(luc[Sample1=="del-uorf_acc",FF_Renilla],luc[Sample1=="uorf_acc",FF_Renilla])
#covUORF<-cov(luc[Sample1=="del-uORF",FF_Renilla],luc[Sample1=="UORF",FF_Renilla])
#covpvORF_acc<-cov(luc[Sample1=="del-pvorf_acc",FF_Renilla],luc[Sample1=="pvorf_acc",FF_Renilla])
#covpvORF<-cov(luc[Sample1=="del-pvORF",FF_Renilla],luc[Sample1=="pvORF",FF_Renilla])

#delEffect<-data.table(comparison=c("uorf_acc","UORF","pvORF","pvorf_acc"),
                      proportion=c(meansLuc[Sample1=="del-uorf_acc",meanFR]/meansLuc[Sample1=="uorf_acc",meanFR],meansLuc[Sample1=="del-uORF",meanFR]/meansLuc[Sample1=="UORF",meanFR]),
                      Error=c(sqrt((meansLuc[Sample1=="del-uorf_acc",sdFR]/meansLuc[Sample1=="del-uorf_acc",meanFR])**2)+((meansLuc[Sample1=="uorf_acc",sdFR]/meansLuc[Sample1=="uorf_acc",meanFR])**2)+2*(covUORF_acc/(meansLuc[Sample1=="del-uorf_acc",meanFR]*meansLuc[Sample1=="uorf_acc",meanFR])),
                              sqrt((meansLuc[Sample1=="del-uORF",sdFR]/meansLuc[Sample1=="del-uORF",meanFR])**2)+((meansLuc[Sample1=="UORF",sdFR]/meansLuc[Sample1=="UORF",meanFR])**2)+2*(covUORF/(meansLuc[Sample1=="del-uORF",meanFR]*meansLuc[Sample1=="UORF",meanFR]))
                      ),
                      proportion2=c(meansLuc[Sample1=="del-pvorf_acc",meanFR]/meansLuc[Sample1=="pvorf_acc",meanFR],meansLuc[Sample1=="del-pvORF",meanFR]/meansLuc[Sample1=="pvORF",meanFR]),
                      Error2=c(sqrt((meansLuc[Sample1=="del-pvorf_acc",sdFR]/meansLuc[Sample1=="del-pvorf_acc",meanFR])**2)+((meansLuc[Sample1=="pvorf_acc",sdFR]/meansLuc[Sample1=="pvorf_acc",meanFR])**2)+2*(covpvORF_acc/(meansLuc[Sample1=="del-pvorf_acc",meanFR]*meansLuc[Sample1=="pvorf_acc",meanFR])),
                               sqrt((meansLuc[Sample1=="del-pvORF",sdFR]/meansLuc[Sample1=="del-pvORF",meanFR])**2)+((meansLuc[Sample1=="pvORF",sdFR]/meansLuc[Sample1=="pvORF",meanFR])**2)+2*(covpvORF/(meansLuc[Sample1=="del-pvORF",meanFR]*meansLuc[Sample1=="pvorf",meanFR]))
                      ))




#ggplot(delEffect)+geom_bar(aes(y=proportion-1,x=comparison,fill=comparison),stat='identity')+geom_errorbar(aes(x=comparison,ymin=proportion-1-Error,ymax=proportion-1+Error),width=0.1)+
  geom_bar(aes(y=proportion2-1,x=comparison2,fill=comparison2),stat='identity')+geom_errorbar(aes(x=comparison2,ymin=proportion2-1-Error2,ymax=proportion2-1+Error2),width=0.1)+
  scale_fill_manual(values = c("cyan", "blue","pink","orange"))+
  theme_pubr()

  