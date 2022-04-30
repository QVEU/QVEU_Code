#!/usr/bin/env bash
touch blasrOutput.out
out=blasrOutput.out
readsDir=/Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fastq_pass/
> $out
for i in $readsDir/*.fastq
do
blasr $i Homo_sapiens.GRCh38.cdna.all.fa >> $out
done
