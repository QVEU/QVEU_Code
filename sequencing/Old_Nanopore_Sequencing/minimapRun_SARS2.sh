#!/usr/bin/env bash

#Run-specific Parameters...
out=minimapOutput_NanoporeRun_00001
target=MT276325.1.fa
readsDir=/Library/MinKNOW/data/SARS2/A549_SCoV2_totRNA/20200731_2123_MN34223_FAO33670_43d615e8/fastq_pass

#Script...
echo "Mapping reads..."
for i in ${readsDir}/*.fastq
do
  echo $i
  /Users/ptdolan/Research\ Projects/SARS-COV2/Nanopore_Sequencing/minimap2/minimap2 -ax splice -uf -k14 $target $i > ${i}.sam  # noisy Nanopore Direct RNA-seq
done

#Post process pile of sams.
echo "Merging and removing sams..."
samtools merge -f ${out}.sam ${readsDir}/*.sam
rm $readsDir/*.sam

echo "Making bam output..."
samtools view -S -b ${out}.sam > ${out}.bam

samtools sort ${out}.bam > ${out}_sort.bam
samtools mpileup ${out}_sort.bam > ${out}.pile
samtools bedcov -Q 8 MT276325.1.bed ${out}.sam
echo "Done. "
