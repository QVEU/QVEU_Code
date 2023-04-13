#!/bin/sh
#$ -l gpu=1,v100

module load guppy_gpu/6.0.6
guppy_basecaller -i /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0057_Insertional_Handle_P1_Replication_Proteins/wbrepproteinsev71handle/no_sample/20230303_2128_MC-113212_FAV68745_76ebd27f//fast5_pass/barcode08 -s /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0057_Insertional_Handle_P1_Replication_Proteins/wbrepproteinsev71handle/no_sample/20230303_2128_MC-113212_FAV68745_76ebd27f/guppy_basecalled/barcode08 --flowcell FLO-MIN106 --kit SQK-NBD112-24 -x cuda:0
