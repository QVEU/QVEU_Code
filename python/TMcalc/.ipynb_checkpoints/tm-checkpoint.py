from Bio.SeqUtils import MeltingTemp
import sys

def Tm(primerSeq):
    TM=MeltingTemp.Tm_NN(primerSeq,Na=50, K=50, Mg=1.5,nn_table=MeltingTemp.DNA_NN4)
    return(TM)

if __name__=="__main__":
    oligo=sys.argv[1]
    print(round(Tm(oligo),2))