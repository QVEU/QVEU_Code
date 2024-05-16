#extractMRCAs.R
#Script to extract mrca table information from BASR csv. 
#Usage: Rscript extractMRCAs.R <BASRcsv.csv>
#

library(data.table)

infile=commandArgs(trailingOnly=TRUE)[1]
lineages=fread(infile)

lineageA<-lineages
lineageB<-lineages

pairwiseMERGE<-merge.data.table(lineageA,lineageB,by=c('V1','i','j'),allow.cartesian = T)
pairwiseMERGE[,lineagePair:=paste(taxLineage.x,taxLineage.y,sep="_")]
pairwiseMERGE[,maxI:=max(i),by=c("lineagePair","V1")]

mrcas=pairwiseMERGE[i==maxI,]

write.csv(mrcas,replace(infile,list = ".csv","_MRCA_Data.csv"))