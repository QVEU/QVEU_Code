#!/bin/sh
#$ -l gpu=1,v100

module load guppy_gpu/6.0.6
guppy_basecaller -i inputdirectory -s savedirectory --flowcell FLO-MIN106 --kit SQK-NBD112-24 -x cuda:0
