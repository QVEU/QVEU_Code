#!/bin/bash

# Job name
#$ -N minimap2 

# Execute the script form the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j n

# Send the output of the script to a directory called 'UGE-output' in the 
# current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE 
# scheduler will NOT create this direcotry for you**
#$ -o /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0027_Miseq_EVD68_Passage5/

# Tell the job your memory requirements
#$ -l mem_free=10G,h_vmem=12G

# Send mail when the job is submitted, and when the job completes
#$ -m be

# Specify an email address to use
#$ -M adam.catching@naiad.nih.gov

#==========================================================================#
# This script is for automating the minimap2 alignment of Illumina reads
# And for batching the reference sequence by row of a 96-well plate


module load minimap2/2.10

# Define base directory
base_dir=/hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0027_Miseq_EVD68_Passage5/

# Define location of reference sequences
fermon_ref="fermon.fa"
MO_ref="MO.fa"

# Loop through N samples
N=96
# for i in {1..24} --> This is old code for reference
for i in $( eval echo {0..$N} )
do
  
  # This function groups the i-th sample into groups of 8 (for reference sequence)
  alternate=$(( ( (i-1)/8 ) % 2))

  # There is an extra 0 in the name of the first 9 reads that will mess up location names
  if (($i < 10))
    then 
      udp="UDP000"$i
    else 
      udp="UDP00"$i
  fi

  # Define the location of the reads for this barcode
  read1=$base_dir"fastq/"$udp"_S"$i"_L001_R1_001.fastq.gz"
  read2=$base_dir"fastq/"$udp"_S"$i"_L001_R2_001.fastq.gz"

  # Define the location of the output files (requires > to be in the string)
  aligned_output=$base_dir"sam/"$udp".sam"

  # Align based on i-th reference sequence
  if (($alternate==0))
      then
        echo minimap2 -ax sr fermon.fa $read1 $read2 ">" $aligned_output
        minimap2 -ax sr $fermon_ref $read1 $read2 > $aligned_output
      else
        echo minimap2 -ax sr MO.fa $read1 $read2 ">" $aligned_output
        minimap2 -ax sr $MO_ref $read1 $read2 > $aligned_output
  fi
done

echo All reads have been mapped