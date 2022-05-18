#Basic Alignment.R
#Alignment of fa files. Usage (in command line): `Rscript BasicAlignment.R <myfasta.fa> [<other optional .fa's ... >]
library(Biostrings)
library(DECIPHER)
inputList=commandArgs(trailingOnly=TRUE)
inputSeqs<-readDNAStringSet(inputList)#Vector of FASTA files.


MED=median(width(inputSeqs))#median length
print(MED)
MAD=mad(width(inputSeqs))#median deviation

prelength<-length(inputSeqs)
inputSeqs<-inputSeqs[abs(width(inputSeqs)-MED)<=(MAD*2)]#filter out long or short sequences (2 deviations away...)
postlength<-length(inputSeqs)

print(paste("Cut sequences from",prelength,"to",postlength,"."))

inputSeqs<-inputSeqs[]#Filter by some criteria. i
names(inputSeqs)<-lapply(X = names(inputSeqs),function(X){paste(strsplit(X," ")[[1]],collapse= "_")}) #avoids some name parsing issues with other programs. 
inputSeqs_Gapless<-RemoveGaps(inputSeqs)
MAlign<-AlignSeqs(processors = 8,inputSeqs_Gapless,refinements=10)
SMAlign<-StaggerAlignment(MAlign)
writeXStringSet(SMAlign,filepath = paste(inputList[1],"_basic.aln",sep = ""))

#system(command = "FastTree -nt -gtr < ./basicAln_output.aln > basicAln.tree")

