library(data.table)
library(ggplot2)
library(ggforce)
# devtools::install_github("stephenturner/annotables")
# library(annotables)

#
Humidor<-fread("/Users/ptdolan/Research Projects/SARS-COV2/Nanopore_Sequencing/Humidor.csv")
Depth<-fread("~/Research_Projects/SARS-COV2/Nanopore_Sequencing/minimapOutput_MT276325.pile",header = F)

#Import Position Annotation
genomeBed<-fread("~/Research_Projects/SARS-COV2/Nanopore_Sequencing/MT276325.1.bed")
genomeBed$index=1:nrow(genomeBed)
genomeBed[,name:=limma::strsplit2(V4,":")[,2]]

#BaseCallFrequencies
CharCounts<-Depth[,data.table(table(toupper(limma::strsplit2(V5,"")))),by=V2]
paste(CharCounts[V1%in%c("A","C","T","G"),V1[N==max(N)],by=V2]$V1,collapse = "")
minimap<-fread(nrows = 10000,"~/Research_Projects/SARS-COV2/Nanopore_Sequencing/minimapOutput_SARS2.sam",comment="@")

#Import Cigar data
Humidor<-fread("/Users/ptdolan/Research_Projects/SARS-COV2/Nanopore_Sequencing/Humidor.csv")

#plot base call freqs
RANGE=1:30000
IDs<-ggplot(CharCounts[(V2%in%RANGE)&(V1%in%c("A","T","C","G","U"))])+geom_bar(aes(V2,N,fill=V1),stat="identity",pos="stack")+xlab("genome position")
ggsave("basecalldepth.pdf",IDs)
RANGE=28200:28950
IDs<-ggplot(CharCounts[(V2%in%RANGE)&(V1%in%c("A","T","C","G","U","<"))])+geom_bar(aes(V2,N,fill=V1),stat="identity",pos="fill")
ggsave("basecallzoom_high.pdf",IDs)

D=ggplot(Humidor[type=="N"])+
  geom_bin2d(aes(refPos,qPos),binwidth=200)+scale_fill_viridis_c(trans="log10")
ggsave("junctionDensity.pdf",D)
  
ggplot(Humidor[type=="N"][1:1000000])+
  geom_point(aes(refPos,qPos),cex=0.5,alpha=0.2)
  
ggplot(Humidor[type=="N"][1:1000000])+  
  geom_histogram(aes(refPos),col="",binwidth = 10)+scale_y_log10()
SCALE=30000
  ggplot(Humidor[type=="N"][1:10000])+
    geom_vline(xintercept = c(genomeBed$V2),col="green")+
    geom_segment(aes(x=refPos,xend=qPos,y=0, yend=SCALE),alpha=0.1)+
    geom_segment(data=genomeBed,aes(x=V2, xend=V3, y=SCALE+(0.02*SCALE*(index)), yend= SCALE+(0.02*SCALE*(index)), col=name),lwd=3)+
    geom_text(data=genomeBed,aes(x=rowMeans(cbind(V2,V3)), y=SCALE+(0.02*SCALE*(index)),label=name))+
    scale_color_discrete()+theme_minimal()

juncCounts=Humidor[type=="N",length(type),by=list(qPos,refPos)]

wholePlot<-ggplot(juncCounts[V1>1])+
  geom_vline(xintercept = c(genomeBed$V2),col="grey")+
  geom_hline(yintercept = c(genomeBed$V2),col="grey")+
  geom_point(aes(refPos,qPos,size=V1, col=V1),pch=1,alpha=0.4)+
  scale_size_continuous("Junc. Count",trans="log10")+
  scale_color_viridis_c("Junc. Count",trans="log10")+
  theme_minimal()

ggsave(plot=wholePlot,filename = "juncCountplot.pdf",width=10, height=8.5)
ggsave(plot=wholePlot+xlim(25000,NA)+ylim(0,10000),filename = "juncCountplot_zoom.pdf",width=10, height=8.5)

ggplot(Humidor[type=="N"])+
  geom_histogram(aes(refPos),stat="density",col='red')+
  geom_histogram(aes(qPos),stat="density",col='blue')+ylab("Site Frequency")+xlab("Donor/Acceptor Site")
  #geom_arc(data=juncCounts,aes(x0=qPos+refPos,lwd=log10(V1),y0=1,r=(abs(refPos-qPos)/2),start=1.5*pi,end=0.5*pi),alpha=0.4,color='black')
