#!/usr/bin/env bash
#/Users/ptdolan/miniconda3/envs/py37/bin/multi_to_single_fast5 --input_path /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/ --save_path /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/single_reads/ --recursive

#tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/single_reads/ --fastq-filenames /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fastq_pass/*.fastq
#Virus Resquiggle
#tombo resquiggle --rna  --processes 8 --num-most-common-errors 5 /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/ /Users/ptdolan/Research\ Projects/SARS-COV2/Nanopore_Sequencing/Templates_Fastas/MT276325.1.fasta

#Human resquiggle
tombo resquiggle --rna  --overwrite --processes 10 --num-most-common-errors 5 /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/single_reads /Users/ptdolan/Research\ Projects/SARS-COV2/Nanopore_Sequencing/Homo_sapiens.GRCh38.cdna.all.fa

tombo detect_modifications alternative_model --alternate-bases 5mC 6mA --fast5-basedirs /Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fast5_pass/ \
   --statistics-file-basename SARS-CoV-2_humanRNAmethylation \
   --processes 10
