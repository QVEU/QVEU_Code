#!/bin/bash/
name=Input
module load minimap2
module load gatk/4.2.6.0-Java-1.8.0_92
minimap2 -ax sr /hpcdata/lvd_qve/QVEU_Code/sequencing/template_fastas/puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta VP4-DMS-P1-RepA_S4_L001_R*_001.fastq.gz  > output.sam
gatk AnalyzeSaturationMutagenesis -I output.sam -R ./puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta --orf 746-7327 -O $name --min-flanking-length 10
