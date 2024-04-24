#AnalyzeTree
library(raster)
#library(sf)
#install.packages(c("sf","units"))
#library(units)
#library(smoothr)
library(data.table)
library(ggplot2)
hdr<-function(X){
  Xs<-X[order(X)]
  out<-data.frame(ymin=NA,y=NA,yin=NA)
  try(out<-data.frame(ymax=Xs[(.25*length(Xs))],y=Xs[(.5*length(Xs))],ymin=Xs[(.75*length(Xs))]))
  return(out)
}

setwd("~/Desktop/Tric_Phylo_2020_Results/")

n<-c("nTail","cTail","NSL","Apical","Equitorial","Intermediate")

plotSelection<-function(titlestring=""){
  title=titlestring
  fn<-paste("/Users/ptdolan/Dropbox/Phylogenetics/AncChargeRecon/ASR_CCT2020_100M_",titlestring,".csv",sep="")
  print(fn)
  lineages<-fread(fn,nThread = 3,header=T)
  
  #How to recover a matrix of shared MRCA?
  # Its the largest Seq[i,j] where j is max shared i,j
  
  lineages[,ID:=replace(list = "\'",values = "",x = limma::strsplit2(taxLineage,"-")[1])[1],taxLineage]
  lineages[,species:=paste(limma::strsplit2(limma::strsplit2(taxLineage,"-")[1],split = "")[,2:3],collapse=""),taxLineage]
  lineages[,prot:=paste(limma::strsplit2(limma::strsplit2(taxLineage,"-")[1],split = "")[,4:length(limma::strsplit2(limma::strsplit2(taxLineage,"-")[1],split = ""))],collapse=""),taxLineage]
  lineages[,maxdepth:=max(i),by=c("taxLineage","V1")]
  lineages[,brDepth:=-(max(brDist)-brDist),by=c("taxLineage","V1")]
  
  P<-acast(lineages,prot~brDist+gravy+aroma+charge+pI)
  print(P)
  #lineages[,brDepth:=max(brDepth)-brDepth,by=c("taxLineage","V1")]
  lineages<-lineages[prot%in%c("CCT4","CCT5","CPN")]
  plot(ggplot(lineages)+ggtitle(titlestring)+scale_fill_brewer(palette = "Paired")+scale_color_brewer(palette = "Paired")+
    stat_summary_bin(aes(brDist,charge,fill=prot),alpha=0.5,binwidth=0.2,geom="ribbon",col='black',lwd=.3,fun.data=hdr))
    #geom_line(cex=0.2,aes(brDist,charge,group=paste(V1,taxLineage),col=prot),alpha=0.05)
  print(lineages)
  plot(ggplot(lineages)+ggtitle(titlestring)+scale_fill_brewer(palette = "Paired")+scale_color_brewer(palette = "Paired")+
         stat_summary_bin(aes(brDist,aroma,fill=prot),,alpha=0.5,col='black',lwd=.3,binwidth=0.2,geom="ribbon",fun.data=hdr))
         #geom_line(cex=0.2,aes(brDist,aroma,group=paste(V1,taxLineage),col=prot),alpha=0.05)
  
  plot(ggplot(lineages)+ggtitle(titlestring)+scale_fill_brewer(palette = "Paired")+scale_color_brewer(palette = "Paired")+
         stat_summary_bin(aes(brDist,pI,fill=prot),,alpha=0.5,binwidth=0.2,geom="ribbon",col='black',lwd=.3,fun.data=hdr))#+
         #geom_line(cex=0.2,aes(brDist,pI,group=paste(V1,taxLineage),col=prot),alpha=0.05)+
         #facet_wrap(~ID))
  
  plot(ggplot(lineages)+ggtitle(titlestring)+scale_fill_brewer(palette = "Paired")+scale_color_brewer(palette = "Paired")+
         stat_summary_bin(aes(brDist,gravy,fill=prot),,alpha=0.5,binwidth=0.2,geom="ribbon",col='black',lwd=.3,fun.data=hdr))#+
         #geom_line(cex=0.2,aes(brDist,gravy,group=paste(V1,taxLineage),col=prot),alpha=0.05)+
         #facet_wrap(~ID))

  # plot(ggplot(lineages)+ggtitle(titlestring)+scale_color_brewer(palette = "Paired")+
  #        stat_summary_bin(aes(brDepth,charge,col=prot),alpha=0.2)+
  #        geom_line(cex=0.2,aes(brDepth,charge,group=paste(V1,taxLineage),col=prot),alpha=0.05))
  # 
  # plot(ggplot(lineages)+ggtitle(titlestring)+scale_color_brewer(palette = "Paired")+
  #        geom_point(cex=0.1,aes(brDepth,charge,col=prot),alpha=0.2)+
  #        geom_line(cex=0.2,aes(brDepth,charge,group=paste(V1,taxLineage),col=prot),alpha=0.05))
  # 
  # plot(ggplot(lineages)+ggtitle(titlestring)+scale_color_brewer(palette = "Paired")+
  #        geom_point(cex=0.1,aes(brDepth,pI,col=prot),alpha=0.2)+
  #        geom_line(cex=0.2,aes(brDepth,pI,group=paste(V1,taxLineage),col=prot),alpha=0.05))
  # 
  # plot(ggplot(lineages)+ggtitle(titlestring)+scale_color_brewer(palette = "Paired")+
  #        geom_point(cex=0.1,aes(brDepth,gravy,col=prot),alpha=0.2)+
  #        geom_line(cex=0.2,aes(brDepth,gravy,group=paste(V1,taxLineage),col=prot),alpha=0.05))
  # 
  # plot(ggplot(lineages)+ggtitle(titlestring)+scale_color_brewer(palette = "Paired")+
  #      geom_point(cex=0.1,aes(brDepth,aroma,col=prot),alpha=0.2)+
  #       geom_line(cex=0.2,aes(brDepth,aroma,group=paste(V1,taxLineage),col=prot),alpha=0.05))
}
pdf("allout.pdf",width=10,height=10)
for (N in n){plotSelection(N)}
dev.off()

dNdS<-fread("/Users/ptdolan/Desktop/Tric_Phylo_2020_Results/TricCPN_codAln_100M.TricCPN_codAln.dNdS.log",nThread = 3,header=T)
mdNdS<-melt(dNdS[state>10E6],id.vars = 1)

#N-tails
nTail=1:51
#C-tails
cTail=628:655
#Nucleotide Sensing Loop
NSL=200:210
#Apical Domain
apical=264:452
#Equitorial Domain
equitorial=c(52:178,491:627)
#Intermediate Domain
intermediate=c(179:263,453:490)

p<-list(nTail,cTail,NSL,apical,equitorial,intermediate)
n<-c("nTail","cTail","NSL","Apical","Equitorial","Intermediate")

mdNdS$pos<-as.integer(limma::strsplit2(mdNdS[,limma::strsplit2(variable,split = "TricCPN_codAln.dNdS\\[")][,2],split = "\\]")[,1])

for (i in 1:length(n)){
  mdNdS[pos%in%p[[i]],region:=as.factor(n[i])]
}

ggplot(mdNdS)+
  stat_summary_bin(aes(pos,value),geom="ribbon",fill='darkgrey',fun.data = hpd95,binwidth = 1)+
  stat_summary_bin(aes(pos,value),lwd=0.3,col="black",geom="line",fun = median,binwidth = 1)+
  stat_summary_bin(cex=0.3,aes(pos,value,col=region),geom="point",fun = median,binwidth = 1)

ggplot(mdNdS)+
  stat_summary_bin(aes(pos,value),geom="ribbon",col='grey',fun.data = hpd95,binwidth = 10)+
  stat_summary_bin(aes(pos,value),lwd=0.3,col="black",geom="line",fun = median,binwidth = 10)+
  stat_summary_bin(cex=0.3,aes(pos,value,col=region),geom="point",fun = median,binwidth = 10)

ggplot(mdNdS)+geom_boxplot(aes(fill=region,region,value))
ggplot(mdNdS)+geom_boxplot(aes(fill=region,region,value))

input<-seqinr::read.fasta("TricCPN_codAln.fa")
alignBlock<-rbind.data.frame(lapply(input,function(X){
  return(as.integer(factor(seqinr::translate(X),levels=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","Y","Z"))))
}))
alignBlockM<- reshape2::melt(data.frame(pos=1:nrow(alignBlock),alignBlock),id.var=1)
ggplot(alignBlockM,aes(pos,variable))+geom_tile(aes(fill=as.character(value)),show.legend = F)+scale_color_viridis_d()+ylab("Sequence")


lapply(unlist(seqinr::read.fasta("TricCPN_codAln.fa")),FUN = function(X){seqinr::translate(X)})

