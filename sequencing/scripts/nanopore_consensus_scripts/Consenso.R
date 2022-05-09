#!/usr/bin/env Rscript
require(data.table)
require(ggplot2)
require(cowplot)

args = commandArgs(trailingOnly=TRUE)
require(Rsamtools)

bam<-Rsamtools::scanBam(paste(inputDir,"merge.bam",sep=""))
unroll_bam<-bam[[1]]
DTbam<-as.data.table(unroll_bam)
DTbam[,readLength:=width(DTbam[,seq])]
DTbam[,seq:=toString(seq)]
DTbam$readSequence
hist(DTbam$readLength)
toString(DTbam$seq)
#print(args)
#inputDir=args[1]
inputDir="~/Research/minion_data/map_cvb/"

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

a=c()
mutRateCutoff=for(i in 1:3000){
  a<-append(a,mad(log10(sample(mutFreqs[Depth>(medianDepth-2*madDepth)&CHAR!=consensus]$freq,replace = T))))
}

10**(log10(median(mutFreqs[Depth>(medianDepth-2*madDepth)&CHAR!=consensus]$freq))+median(a))

freqPlot=ggplot(mutFreqs[CHAR!=consensus,])+
  geom_hline(yintercept = median(a),lty=2,alpha=0.5)+
  geom_text(x=0,y = median(a),label=round(median(a),2),hjust=T)+
  ggtitle(paste("Reads mapped to: ",inputDir))+
  scale_color_brewer(palette = "Set1")+
  geom_point(show.legend = F,aes(POS,freq,col=Depth>(medianDepth-2*madDepth),alpha=(Depth>(medianDepth-2*madDepth))&(freq>median(a))))+
  geom_text(cex=2,hjust = T ,data = mutFreqs[CHAR!=consensus&(Depth>(medianDepth-2*madDepth))&(freq>median(a)),],aes(POS,freq,label=paste(consensus,POS,CHAR,sep="")))
10**.57

ggsave2(cowplot::plot_grid(freqPlot,covplot,nrow = 2),filename = paste(inputDir,"coverage_and_variantFreq.pdf",sep=""),width = 7,height=7)

#write table of significant variants. 
write.csv(file = paste(inputDir,"variant_table.csv",sep=""),mutFreqs[(CHAR!=consensus)&(Depth>(medianDepth-2*madDepth))&(freq>median(a))])

