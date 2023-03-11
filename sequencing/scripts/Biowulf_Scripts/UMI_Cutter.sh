#!/bin/bash/
#This bash script is used to remove UMIs from NGS data using Twist's UMI.
#It should be run on the biowulf servers.
#1 Picard: FastqToSam
module load picard/2.27.3
java -Xmx2048m -jar $PICARDJARPATH/picard.jar FastqToSam \
       F1=VP4-DMS-Input_S1_L001_R1_001.fastq.gz \
       F2=VP4-DMS-Input_S1_L001_R2_001.fastq.gz \
       O=VP4-DMS-Input_S1_L001_unaligned.bam \
       SM=sample001 \
       RG=rg001
module unload picard/2.27.3

#2 fgbio: ExtractUmisFromBam (5M2ST+)
#Read structure: (5 bases of molecular index then remaining bases) x2
module load fgbio/2.0.2
java -Xmx2048m -jar $FGBIOJAR ExtractUmisFromBam \
-i VP4-DMS-Input_S1_L001_unaligned.bam \
-o VP4-DMS-Input_S1_L001_unaligned_UMIextracted.bam \
-r 5M2S+T 5M2S+T \
--molecular-index-tags RX
module unload fgbio/2.0.2

#3 Picard: SamToFastq
module load picard/2.27.3
java -Xmx2048m -jar $PICARDJARPATH/picard.jar SamToFastq \
            -I VP4-DMS-Input_S1_L001_unaligned_UMIextracted.bam \
            -F VP4-DMS-Input_S1_L001_unaligned_UMIextracted_R1.fastq \
            -F2 VP4-DMS-Input_S1_L001_unaligned_UMIextracted_R2.fastq
module unload picard/2.27.3
