library(data.table)
library(ggplot2)
devtools::install_github("stephenturner/annotables")
library(annotables)
Blasr<-fread("~/Research Projects/SARS-COV2/Nanopore_Sequencing/blasrOutput.out")
colnames(Blasr)<-c("qName","tName","qStrand","tStrand","score","percentSimilarity","tStart","tEnd","tLength","qStart","qEnd","qLength","nCells")
Blasr[,nhits_transcript:=length(qName),by=tName]
Blasr[,gene:=limma::strsplit2(tName,"\\.")[,1]]
Blasr[,nhits_gene:=length(qName),by=gene]
length(Blasr$tName)
length(unique(Blasr$gene))

ENSmerge=merge.data.table(Blasr,annotables::grch38_tx2gene,by.x = "gene",by.y="enstxp")

fullmerge=merge.data.table(ENSmerge,annotables::grch38,by="ensgene")

fullmerge[,tophit:=.SD[.SD$score==min(.SD$score),symbol][1],by=qName]

ggplot(fullmerge[symbol==tophit])+geom_point(aes(start,nhits_gene,col=biotype))+facet_wrap(~chr)+scale_y_log10()

Blasr<-fread("~/Research Projects/SARS-COV2/Nanopore_Sequencing/blasrOutput_SARS2.out")
colnames(Blasr)<-c("qName","tName","qStrand","tStrand","score","percentSimilarity","tStart","tEnd","tLength","qStart","qEnd","qLength","nCells")
Blasr[,nhits_transcript:=length(qName),by=tName]
length(Blasr$tName)
length(unique(Blasr$tName))

Blasr

ggplot(Blasr[1:10000])+
  geom_errorbarh(aes(col=tName,xmin=tStart,xmax=tEnd,y=reorder(qName,tLength)))

minimap<-fread("~/Research Projects/SARS-COV2/Nanopore_Sequencing/minimapOutput_SARS2.sam")
