
library(ggpubr)
luc<-fread("~/Documents/luciferase data.csv")
luc<-fread("~/Documents/RW_luc_data_2.csv")
luc<-fread("~/Documents/RRL_luc2_data.csv")
luc<-fread("~/Documents/RRL_luc_data.csv")
luc=luc[,1:7]#remove analysis columns

luc[,meanFR:=mean(FF_Renilla),by=list(Sample1)]

luc

luc[,sdFR:=sd(FF_Renilla), by=list(Sample1)]

luc

#T-tests
t.test(luc[Sample1=="UORF",FF_Renilla],luc[Sample1=="del-uORF",FF_Renilla])
t.test(luc[Sample1=="pvORF",FF_Renilla],luc[Sample1=="del-pvORF",FF_Renilla])


ggplot(luc[Sample1!="GFP"&Rep=="1A"]) + 
  geom_bar(aes(y=meanFR,x = reorder(Sample1,meanFR,mean),fill=Sample1),stat = "identity")+
  geom_errorbar(aes(ymin=meanFR-sdFR, ymax=meanFR+sdFR,x = Sample1))+
  scale_fill_manual(values = c("cyan", "lightblue", "grey", "blue", "purple"))+
  theme_pubr()+xlab("Sample")+ylab("Firefly/Renilla")



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
                    
  