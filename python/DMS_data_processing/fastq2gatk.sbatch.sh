#!/bin/bash

# This script takes data exported from the MiSeq and converts to an array of DENV2 amino acid frequencies

# Job Name
#$ -N 

#SBATCH --job-name=fastq2gatk
#SBATCH --output=fastq2gatk.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64g

#  Specify an email address to use
#$ -M catchingba@nih.gov

#module load minimap2/2.26-hnckepf
#module load gatk/4.2.6.0-Java-1.8.0_92
module load bcl2fastq2/2.20.0.422-orocbiu


# Make directories of different data formats
mkdir fastq
mkdir sam
mkdir gatk

# Convert .bcl files to .fastq files and barcode from samplesheet
#bcl2fastq -o fastq
#echo Done converting to fastq files
module purge
module load minimap2/2.26-hnckepf
module load gatk/4.5.0.0 
# Define the number of the last barcode 
N=68
# Define the number of the first barcode
n=63
for i in $( eval echo {$n..$N} )
do
  # Define the barcode index (from samplesheet)
  udp="UDP_00"$i
  # ONLY ADD FOR NEXTSEQ RUNS
  s="$(($i-$n+1))"

  # Define the location of the output files (requires > to be in the string)
  aligned_output="sam/"$udp".sam"
  # Define the gatk output
  gatk_output="gatk/"$udp

  # Define the location of the reads for this barcode
  read1="fastq/"$udp"_S"$s"_R1_001.fastq.gz"
  read2="fastq/"$udp"_S"$s"_R2_001.fastq.gz"

  # Align to reference sequence
  minimap2 -ax sr  pUC19_EVD68_US_14_IL_18952.fasta $read1 $read2 > $aligned_output
  echo "aligning barcode"$i

  # Map codons 
  gatk AnalyzeSaturationMutagenesis -I $aligned_output -R  pUC19_EVD68_US_14_IL_18952.fasta --orf 698-7264 -O $gatk_output
  echo "mapping codons of barcode "$i
done
echo Done