#!/usr/bin/env Rscript
require(data.table)
require(ggplot2)
require(cowplot)

args = commandArgs(trailingOnly=TRUE)


inputDir=args[1]
#inputDir="~/Research/minion_data/map_cvb/"

Depth<-fread(paste(inputDir,"merge_sort_pile.pile",sep=""),header = F)

#BaseCallFrequencies
CharCounts<-Depth[,data.table(table(toupper(limma::strsplit2(V5,"")))),by=V2]
colnames(CharCounts)<-c("POS","CHAR","DEPTH")
#write consensus sequence

mutFreqs<-CharCounts[CHAR%in%c("A","C","G","T"),list(CHAR,freq=DEPTH/sum(DEPTH),Depth=sum(DEPTH),consensus=CHAR[DEPTH==max(DEPTH)]),by=POS]

freqTable<-dcast.data.table(mutFreqs,POS+consensus+Depth~CHAR,value.var = "freq",fill = 0.0)
medianDepth<-median(freqTable$Depth)
madDepth<-mad(freqTable$Depth)
CONSENSUS<- freqTable[,ifelse(Depth>(medianDepth-2*madDepth),consensus,"-"),by=POS]

write(file=paste(inputDir,"consenso_consensus.fa",sep=""),x = paste(">consensus\n",paste(CONSENSUS$V1,collapse = "")))

#writeFreqTable
write.csv(file = paste(inputDir,"freqtable.csv",sep = ""),x = freqTable)

covplot = ggplot(freqTable)+scale_color_brewer(palette = "Set1")+
  geom_smooth(stat=mean_cl_boot())+
  geom_line(show.legend = F,aes(POS,Depth,group=NA,col=Depth>(medianDepth-2*madDepth)))

mutMAD=c()
mutRateCutoff=for(i in 1:3000){
  mutMAD<-append(a,mad(sample(mutFreqs[Depth>(medianDepth-2*madDepth)&CHAR!=consensus]$freq,replace = T)))
}

freqPlot=ggplot(mutFreqs[CHAR!=consensus,])+
  geom_hline(yintercept = median(a),lty=2,alpha=0.5)+
  geom_text(x=0,y = median(a),label=round(median(a),2),hjust=T)+
  ggtitle(paste("Reads mapped to: ",inputDir))+
  scale_color_brewer(palette = "Set1")+
  geom_point(show.legend = F,aes(POS,freq,col=Depth>(medianDepth-2*madDepth),alpha=(Depth>(medianDepth-2*madDepth))&(freq>median(mutMAD))))+
  geom_text(cex=2,hjust = T ,data = mutFreqs[(CHAR!=consensus)&(Depth>(medianDepth-2*madDepth))&(freq>median(mutMAD)),],aes(POS,freq,label=paste(consensus,POS,CHAR,sep="")))

ggsave2(plot_grid(freqPlot,covplot,nrow = 2),filename = paste(inputDir,"coverage_and_variantFreq.pdf",sep=""),width = 7,height=7)

#write table of significant variants. 
write.csv(file = paste(inputDir,"variant_table.csv",sep=""),mutFreqs[(CHAR!=consensus)&(Depth>(medianDepth-2*madDepth))&(freq>median(a))])

