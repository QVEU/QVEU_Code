#!/bin/bash

# Job name
#$ -N Sam2Q20

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
# This script is for automating the generation of a Q20-cutoff of aligned EV-D68 genomes

module load python/3.7.3
module load numpy/1.14.0-goolf-1.7.20-Python-2.7.11 

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
      udp="barcode0"$i
    else 
      udp="barcode"$i
  fi

  # Define the location of the reads for this barcode
  samfile="sam/"$udp".sam"

  # Define the location of the output files (requires > to be in the string)
  q20_output="Q20/"$udp".txt"

  # Align based on i-th reference sequence
  if (($alternate==0))
      then
        echo "sam_to_Q20.py"$MO_ref$q20_output
        python sam_to_Q20.py $MO_ref $samfile $q20_output
      else
        echo "sam_to_Q20.py"$fermon_ref$q20_output
        python sam_to_Q20.py $fermon_ref $samfile $q20_output
  fi
done

echo All reads have been mapped