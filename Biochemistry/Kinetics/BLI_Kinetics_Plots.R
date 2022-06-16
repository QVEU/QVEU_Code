library(ggplot2)
ptbfrags<-read.csv("~/GitHub/QVEU_Code/Biochemistry/Kinetics/RW_BLI_2022_06_15_wt_del_rep_ptbfrags.csv")

ggplot(ptbfrags)+
  geom_boxplot(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(KD..M.),fill=`Loading.Sample.ID`))+
  geom_point(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(KD..M.)))+
  ylab(expression(K[D]))+xlab("IRES - PTB Combination")+  
  scale_y_log10()+theme_bw()

ggplot(ptbfrags)+
  geom_boxplot(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(ka..1.Ms.),fill=`Loading.Sample.ID`))+
  geom_point(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(ka..1.Ms.),fill=`Loading.Sample.ID`))+
  ylab(expression(K[a]))+xlab("IRES - PTB Combination")+  
  scale_y_log10()+theme_bw()

ggplot(ptbfrags)+
  geom_boxplot(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(kdis..1.s.),fill=`Loading.Sample.ID`))+
  geom_point(aes(paste(Experiment,Loading.Sample.ID, sep = " + "),as.numeric(kdis..1.s.)))+
  ylab(expression(K[dis]))+xlab("IRES - PTB Combination")+  
  scale_y_log10()+theme_bw()

ggplot(ptbfrags)+
  #geom_segment(aes(y=kdis..1.s., yend=kdis..1.s.,x=(ka..1.Ms.-ka.Error), xend=(ka..1.Ms.+ka.Error)))+
  #geom_segment(aes(y=kdis..1.s.-kdis.Error, yend=kdis..1.s.+kdis.Error,x=(ka..1.Ms.), xend=(ka..1.Ms.)))
  geom_point(aes(as.numeric(ka..1.Ms.),as.numeric(kdis..1.s.),color=Loading.Sample.ID,shape=Experiment,size=log10(Conc...nM.)))+
  labs(colour="IRES - PTB Combination",shape="IRES Variant")+
  ylab(expression(K[dis]))+xlab(expression(K[a]))+  
  scale_y_log10()+scale_x_log10()+theme_bw()

