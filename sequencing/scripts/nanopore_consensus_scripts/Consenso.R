library(data.table)
library(ggplot2)
#library(ggforce)
# devtools::install_github("stephenturner/annotables")
# library(annotables)

#CHANGE THIS to your output directory from the mapping script (where your pile up lives). 
#inputDir="~/Research/minion_data/map_evd68/"
inputDir="~/Research/minion_data/map_cvb/"

Depth<-fread(paste(inputDir,"merge_sort_pile.pile",sep=""),header = F)

#BaseCallFrequencies
CharCounts<-Depth[,data.table(table(toupper(limma::strsplit2(V5,"")))),by=V2]
CharCounts
colnames(CharCounts)<-c("POS","CHAR","DEPTH")
#write consensus sequence
write(file=paste(inputDir,"consenso_consensus.fa",sep=""),x = paste(">consensus\n",paste(CharCounts[CHAR%in%c("A","C","T","G"),CHAR[DEPTH==max(DEPTH)],by=POS]$V1,collapse = "")))

freqTable=dcast.data.table(CharCounts[CHAR%in%c("A","C","G","T"),list(CHAR,freq=DEPTH/sum(DEPTH),Depth=sum(DEPTH),consensus=CHAR[DEPTH==max(DEPTH)]),by=POS],POS+consensus+Depth~CHAR,value.var = "freq",fill = 0.0)
#writeFreqTable
write.csv(file = paste(inputDir,"freqtable.csv",sep = ""),x = freqTable)

ggsave(paste(inputDir,"coveragePlot.pdf",sep = ""),plot = ggplot(freqTable)+geom_line(aes(POS,Depth)))

