#!/usr/bin/env bash
#minimap_nanopore_consensus.sh
#usage: sh minimap_nanopore_consensus.sh templatefasta.fa input.fastq
/Users/dolanpt/opt/anaconda3/pkgs/minimap2-2.24-h1f540d2_1/bin/minimap2 -ax map-ont $1 $2 > ${2/\.fastq/}.sam      # for Oxford Nanopore reads
