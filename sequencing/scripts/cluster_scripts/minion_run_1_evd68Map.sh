#!/bin/bash

#nanopore_consensus_run_<yourinfo>.sh
#edit this script for individual runs on the cluster.
#open an interactive job with 'qrsh' and then sun `sh minion_run_<yourinfo>.sh`.

#load required modules...
module load minimap2
module load samtools

#run minimap and samtools to generate consensus
sh ~/QVEU_Code/sequencing/scripts/nanopore_consensus_scripts/minimap_minion_consensus.sh ~/QVEU_Code/sequencing/template_fastas/puc19-evd68-fermon-inverted-2589-9956.fasta ~/Nanopore_Data/fastq_pass/fastqs/
