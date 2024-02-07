#!/bin/bash

# This script performs the Kallisto pipeline on the input fastq files

# Job Name
#$ -N kallisto

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
# **NOTE: Be sure to create this directory before submitting your job, UGE scheduler will NOT create this**
# **directory for you**
#$ -o /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0109_NextSeq_BAC_Fermon_RD_A549/231208_VH01023_36_AAC2LGHM5

# Tell the job your memory requirements
#$ -l mem_free=24G,h_vmem=72G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M YOUR_ADDRESS_HERE

# Load module
module load kallisto/0.45.0-goolf-1.7.20

# Create index of your reference genome (Do this only once), index fasta files can be found here: https://feb2023.archive.ensembl.org/index.html
# kallisto index -i human Homo_sapiens.GRCh38.cdna.all.fa.gz

# Run kallisto transcript quantification, using 30 bootstraps, with the reference human index computed above, to the output folder: WT_RD_37
# Make sure to pair the input fastq.gz files with R1/R2 of each replicate

kallisto quant -b 30 -i human -o WT_RD_37 WT_RD_37_1_S13_R1_001.fastq.gz WT_RD_37_1_S13_R2_001.fastq.gz
