#Basic Alignment.R
#Alignment of fa files. Usage (in command line): `Rscript BasicAlignment.R <myfasta.fa> [<other optional .fa's ... >]

suppressMessages(library(Biostrings))
suppressMessages(library(DECIPHER))

inputList<-commandArgs(trailingOnly=TRUE)
inputSeqs<-readAAStringSet(inputList)#Vector of FASTA files.

inputSeqs<-RemoveGaps(inputSeqs)#Dealign any alignments for correct sequence length.

MED=median(width(inputSeqs))#median length
MAD=mad(width(inputSeqs))#median deviation

prelength<-length(inputSeqs)
#inputSeqs<-inputSeqs[abs(width(inputSeqs)-MED)<=(MAD*2)]#filter out long or short sequences (2 deviations away...)
postlength<-length(inputSeqs)

print(paste("Cut sequences from",prelength,"to",postlength,"."))

inputSeqs<-inputSeqs[]#Filter by some criteria.
names(inputSeqs)<-lapply(X = names(inputSeqs),function(X){paste(strsplit(X," ")[[1]],collapse= "_")}) #avoids some name parsing issues with other programs.

MAlign<-AlignSeqs(processors = 8,inputSeqs,refinements=10)
SMAlign<-StaggerAlignment(MAlign)
writeXStringSet(SMAlign,filepath = paste(inputList[1],"_basicAA.aln",sep = ""))
print(paste("Finished Alignment. Wrote",inputList[1],"_basicAA.aln",sep=""))
