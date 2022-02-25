#delPrime.py: script for designing deletional primer pools
#author: Patrick T. Dolan, QVEU-NIAID-NIH
#USAGE: python3 delPrime.py ./sequenceFileAsSingleLine.txt OutputFileName.txt <deletionSize int (e.g. 3)> <startSite (1 indexed)> <endSite (1 indexed)> 
#Imports

from Bio.SeqUtils import MeltingTemp
from Bio.Seq import Seq
import sys

#Functions

''' Function: Tm()
        Takes primer sequence and compute melting temp. 
    Args:     primerSeq (string) 
    Returns:  melting temperature (float)
'''
def Tm(primerSeq):
    TM=MeltingTemp.Tm_NN(primerSeq,Na=50, K=50, Mg=1.5,nn_table=MeltingTemp.DNA_NN4)
    return(TM)

''' Function: seqWalk()
        Takes primer sequence and compute melting temp. 
    Args:     
        inputSeq: sequence for generating deletions
        delLength: length of deletion (for example,  3 for single codons)
        start: first nt position of first deletion
        end: first nt position of last deletion
    
    Returns:  all the positions of the deletions
'''
def seqWalk(inputSeq,delLength,start=0,end=0):
    if end==0:
        end=len(inputSeq)-1
    positions=range(int(len(inputSeq)/delLength))
    positionsIter=[i*delLength for i in positions if (i >= start/3 and i <= end/3)]
    return(positionsIter)

''' Function: windowSizing()
        Takes primer sequence and compute melting temp. 
    Args:     
        position: site of deletion opening
        delLength: length of deletion (for example,  3 for single codons)
        inputSeq: sequence for generating deletions
        tmDes: desired Tm (minimum)
        windowSize: minimum primer length (will add bases until TmDes is achieved.
    Returns:  outputString - string for output csv file. 
'''
def windowSizing(position,delLength,inputSeq,tmDes,windowSizeMin):
    #reversePrimer (left)
    inputSeq= inputSeq.upper()
    TM=0.0
    windowSize=windowSizeMin
    while TM<tmDes:
        primerR_for=inputSeq[(position-windowSize-2):position-2]
        TM=Tm(primerR_for)
        windowSize+=1
        forSeq = Seq(primerR_for)
        primerR = str(forSeq.reverse_complement())
        
    TM_r=TM #final Tm of reverse annealing region.
    
    #forwardPrimer (right)
    TM=0.0
    windowSize=windowSizeMin
    while TM<=tmDes:
        primerF=inputSeq[(position+delLength-2):(windowSize+position-2)]
        TM=Tm(primerF)
        windowSize+=1
    
    TM_f=TM #final Tm of forward annealing region.
    
    #overlap
    TM=0.0
    windowSize=windowSizeMin
    overlap=5
    while TM<(tmDes+3):
        Fflank=primerR_for[(len(primerR_for)-overlap):len(primerR_for)]
        delPrimerF=(Fflank.lower()+primerF)
        
        Rflank=primerF[0:overlap]
        delPrimerR_for=primerR_for+Rflank.lower()
        
        forSeq = Seq(delPrimerR_for)
        delPrimerR = str(forSeq.reverse_complement())
        overlap+=1
        
        TM=Tm(delPrimerF[0:2*(overlap)])
        
    TM_OL=TM #final Tm of overlap.
    
    print("Deletion at pos: "+str(position)+" to "+str(position+delLength-1))
    print("Forward Primer: "+delPrimerF+"\tLength: "+str(len(primerF))+"\tTm: "+str(TM_f))
    print("Reverse Primer: "+delPrimerR+"\tLength: "+str(len(primerR))+"\tTm: "+str(TM_r))
    print("TM_OL: "+str(TM_OL))
    
    outputString="del_"+str(position)+"-"+str(position+delLength-1)+"_F,"+delPrimerF+","+str(round(TM_f,2))+"\n"+"del_"+str(position)+"-"+str(position+delLength-1)+"_R, "+delPrimerR+","+str(round(TM_r,2))+"\n"#string for text file (e.g. benchling input)
    return(outputString)



#------------------------------------------------------------------
#-------------------------------MAIN-------------------------------
#------------------------------------------------------------------

if __name__ == "__main__":
    #defaults
    tmDes=56
    windowSize=8

    seqFile=sys.argv[1]
    OF=sys.argv[2]
    delLength=int(sys.argv[3])
    start=int(sys.argv[4])
    end=int(sys.argv[5])
    
    with open(seqFile, "r") as IF:
        inputSequence="".join([line for line in IF.readlines()])

    with open(OF,"w") as outfile:
        posArray=seqWalk(inputSequence,delLength,start=start,end=end)
        outfile.write("name,sequence,Tm\n")#header
        for i in posArray:
            outputString=windowSizing(i, delLength, inputSequence, tmDes, windowSize)
            outfile.write(outputString)
            