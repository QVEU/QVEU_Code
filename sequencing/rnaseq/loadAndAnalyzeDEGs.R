library(data.table)
library(stringr)
library(ggplot2)
library(reshape2)
listOfFiles<-list.files("/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0005_NovaSeq_RD+3EVs/",recursive = T,pattern = "differential_expression.csv")

setwd("/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0005_NovaSeq_RD+3EVs/")
allDTs<-data.table()
for (i in listOfFiles[endsWith(listOfFiles,"graphclust/differential_expression.csv")]){
  #print(i)
  name=strsplit(i,"/")[[1]][1]
  print(name)
  thing<-fread(i)
  thing[,sampleID:=name]
  allDTs<-rbindlist(list(allDTs,thing),fill = TRUE)
}
allDTs
mDT<-na.omit(melt.data.table(allDTs,id.vars = c("Feature ID","Feature Name","sampleID"),variable.factor = T))
mDT[,cluster:=str_split(variable," ")[[1]][2],by="variable"]
mDT[,stat:=paste(str_split(variable," ")[[1]][c(3)],str_split(variable," ")[[1]][c(4)],sep = "_"),by="variable"]

ClusterFeatures<-dcast.data.table(mDT,formula = `Feature ID`+`Feature Name`+cluster+sampleID~stat,value.var = "value")

ClusterFeatures[Adjusted_p<0.0000001,Adjusted_p:=0.0000001]

ClusterFeatures[`Feature ID`%in%c("CVB3_Nancy","EVA71-Tainan","EVD68-MO"),rowSums(.SD),by=cluster]

ggplot(ClusterFeatures[`Feature ID`%in%c("CVB3_Nancy","EVA71-Tainan","EVD68-MO")])+
  geom_point(aes(x = cluster,y=Mean_Counts,cex=-log10(Adjusted_p),color=Adjusted_p<0.01))+
  facet_wrap(~sampleID)

ggplot(ClusterFeatures[`Feature ID`%in%c("CVB3_Nancy","EVA71-Tainan","EVD68-MO")])+
  geom_point(aes(x = cluster,y=Log2_fold,cex=-log10(Adjusted_p),color=Adjusted_p<0.01))+
  facet_wrap(~sampleID)

ggplot(ClusterFeatures[Adjusted_p<0.1])+
  geom_hline(yintercept = -log10(0.01))+
  geom_point(cex=0.2,show.legend = F,aes(Log2_fold,-log10(Adjusted_p),color=ClusterFeatures[Adjusted_p<0.1]$`Feature ID`%in%c("CVB3_Nancy","EVA71-Tainan","EVD68-MO")))+
  facet_grid(as.integer(cluster)~sampleID)

DC<-acast(ClusterFeatures,`Feature ID`~cluster+sampleID,value.var = "Mean_Counts")

DCnao<-na.omit(DC)

heatmap(cor(DCnao),scale = "none")
DC[is.na(DC)]<-0

pcDC<-prcomp(DC)
distDC<-cmdscale(dist(t(DC),method = "manhattan"))

#ggplot(data.frame(pcDC$rotation))+geom_point(aes(PC1,PC2))

ggplot(data.frame(distDC))+geom_point(aes(X1,X2,color=log10(DC["CVB3_Nancy",])))
ggplot(data.frame(distDC))+geom_point(aes(X1,X2,color=log10(DC["EVD68-MO",])))
ggplot(data.frame(distDC))+geom_point(aes(X1,X2,color=log10(DC["EVA71-Tainan",])))

ggplot(ClusterFeatures)