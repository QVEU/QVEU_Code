#inPrime.py: script for designing insertional primer pools
#author: Patrick T. Dolan, QVEU-NIAID-NIH
#USAGE: python3 inPrime.py ./sequenceFileAsSingleLine.txt OutputFileName.txt <insertion nt sequence> <startSite (1 indexed)> <endSite (1 indexed)> 

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
def seqWalk(inputSeq,interval=3,start=0,end=0):
    if end==0:
        end=len(inputSeq)-1
    positions=range(int(len(inputSeq)/interval))
    positionsIter=[i*interval for i in positions if (i >= start/interval and i <= end/interval)]
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
def windowSizing(position,interval,inputSeq,insertionSeq,tmDes,windowSizeMin):
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
        primerF=inputSeq[(position-2):(windowSize+position-2)]
        TM=Tm(primerF)
        windowSize+=1
    
    TM_f=TM #final Tm of forward annealing region.
    
    #linker annealing
    linkerSeq_F=primerR_for.lower()+insertionSeq.upper()+primerF.lower()
    linkerSeq_R=str(Seq(linkerSeq_F).reverse_complement())
    
    print("Insertion at pos: "+str(position))
    print("Forward Primer: "+primerF+"\tLength: "+str(len(primerF))+"\tTm: "+str(TM_f))
    print("Reverse Primer: "+primerR+"\tLength: "+str(len(primerR))+"\tTm: "+str(TM_r))

    print("Forward Linker: "+linkerSeq_F+"\tLength: "+str(len(linkerSeq_F)))
    print("Reverse Linker: "+linkerSeq_R+"\tLength: "+str(len(linkerSeq_R)))
    outputString="ins_"+str(position)+"_F,"+primerF+","+str(round(TM_f,2))+"\n"+"ins_"+str(position)+"_R, "+primerR+","+str(round(TM_r,2))+"\n"+"ins_"+str(position)+"_linker_F, "+linkerSeq_F+"\n"+"ins_"+str(position)+"_linker_R, "+linkerSeq_R+"\n"#string for text file (e.g. benchling input)
    return(outputString)

#------------------------------------------------------------------
#-------------------------------MAIN-------------------------------
#------------------------------------------------------------------

if __name__ == "__main__":
    #defaults
    tmDes=56
    windowSize=20 #minimum priming region length (must be 20 for TWIST...)
    
    seqFile=sys.argv[1]
    OF=sys.argv[2]
    interval=int(sys.argv[3])
    
    FLAG="GACTACAAAGACGATGACGACAAG"     #DYKDDDDK
    FLAGx3="GACTATAAAGACCATGATGGAGATTATAAGGACCATGACATCGACTATAAGGATGACGACGATAAA"#DYKDHDG-DYKDHDI-DYKDDDDK
    SIINFEKL="TCCATTATTAATTTCGAGAAACTG" #SIINFEKL
    SUN_GCN4="tcgggcgaggagttattgagcaaaaactatcatttagagaacgaagtcgcgcgcttaaagaaaggctcgggc" #SG-EELLSKNYHLENEVARLKK-GSG
    SUN_GCNx5 = "tcgggtgaagaattactgagtaaaaattatcatctggaaaatgaggtagcgagactaaaaaaggggagtggttctggcgaagagttgctatcgaaaaattatcatcttgagaacgaagttgctaggctcaaaaagggctcaggctcaggcgaggagttgctctcgaaaaactaccacttggaaaatgaggtcgcgaggttgaaaaaggggagcgggtcgggcgaggagttattgagcaaaaactatcatttagagaacgaagtcgcgcgcttaaagaaaggctcgggctcgggcgaagaactcttatcgaagaactaccacctcgaaaatgaggtcgccaggttgaaaaagggcagtggc" # SG-EELLSKNYHLENEVARLKK-GSG x 5
    
    insertionSeq=SIINFEKL
    start=int(sys.argv[4])
    end=int(sys.argv[5])
    
    with open(seqFile, "r") as IF:
        inputSequence="".join([line for line in IF.readlines()])

    with open(OF,"w") as outfile:
        posArray=seqWalk(inputSequence,interval,start=start,end=end)
        outfile.write("name,sequence,Tm\n")#header
        for i in posArray:
            outputString=windowSizing(i, interval, inputSequence, insertionSeq, tmDes, windowSize)
            outfile.write(outputString)