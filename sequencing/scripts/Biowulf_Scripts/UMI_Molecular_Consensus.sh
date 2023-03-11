#!/bin/bash/
#This bash script is used to analyze NGS data using Twist's universal
#adapter system with UMIs The bashscript will take the two fastq files
#originating from sequencing as the initial input. It will extract the
# UMI information from sequencing reads and collapse reads having the UMIs
#It should be run on the biowulf servers
#
#
#
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

#4 BWA: bwa mem (Align reads to reference genome)
module load bwa-mem2/2.2.1
#index your reference FASTA (important step)
bwa-mem2 index puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta
#Align your fastq/-t is for number of threads: 10 set as default
bwa-mem2 mem -t 10 puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta \
VP4-DMS-Input_S1_L001_unaligned_UMIextracted_R1.fastq  VP4-DMS-Input_S1_L001_unaligned_UMIextracted_R2.fastq \
> VP4-DMS-Input_S1_L001_aligned_UMIextracted.bam
module unload bwa-mem2/2.2.1

#5 Merge aligned sequences and unaligned BAM tag files
module load picard/2.27.3
java -Xmx2048m -jar $PICARDJARPATH/picard.jar MergeBamAlignment \ #being killed check memory?
-ALIGNED VP4-DMS-Input_S1_L001_aligned_UMIextracted.bam \
-UNMAPPED VP4-DMS-Input_S1_L001_unaligned_UMIextracted.bam \
-O merge_alignments_VP4-DMS-Input.bam \
-R puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta

#6 Mark duplicates for un-collapsed reads (Quality check step only)
java -Xmx2048m -jar $PICARDJARPATH/picard.jar MarkDuplicates \
 I=merge_alignments_VP4-DMS-Input.bam \
 O=merge_alignments_marked_duplicates_VP4-DMS-Input.bam \
 M=marked_dup_metrics.txt

#7 Group reads by UMI
module load fgbio/2.0.2
java -Xmx2048m -jar $FGBIOJAR GroupReadsByUmi \
-i merge_alignments_VP4-DMS-Input.bam \
-o umi_grouped_VP4-DMS-Input.bam \
-s paired \
-t RX \
-f family_size_counts.txt

#8 fgbio: CallDuplexConsensusReads

java -Xmx2048m -jar $FGBIOJAR CallMolecularConsensusReads \
-i umi_grouped_VP4-DMS-Input.bam \
-o called_umi_grouped_VP4-DMS-Input.bam \
--min-reads 2

#9

module load picard/2.27.3
java -Xmx2048m -jar $PICARDJARPATH/picard.jar SamToFastq \
            -I called_umi_grouped_VP4-DMS-Input.bam \
            -F umi_grouped_VP4-DMS-Input_R1.fastq \
            -F2 umi_grouped_VP4-DMS-Input_R2.fastq
module unload picard/2.27.3

module load bwa-mem2/2.2.1
bwa-mem2 mem -t 10 puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta \
umi_grouped_VP4-DMS-Input_R1.fastq  umi_grouped_VP4-DMS-Input_R2.fastq \
> UMIgrouped_VP4-DMS-Input_Final.bam
module unload bwa-mem2/2.2.1
