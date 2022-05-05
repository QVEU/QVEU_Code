#!/usr/bin/env bash

#Run-specific Parameters...
out=minimapOutput_FINDCOVID
target=/Users/ptdolan/Documents/Research_Projects/SARS-COV2/NC_045512v2.fa
readsDir=/Users/ptdolan/Downloads/FIND_COVID/

#Script...
echo "Mapping reads..."
for i in ${readsDir}/*.fastq
do
  echo $i
  /Users/ptdolan/Documents/Research_Projects/SARS-COV2/Nanopore_Sequencing/minimap2/minimap2 -ax splice -uf -k14 $target $i > ${i}.sam  # noisy Nanopore Direct RNA-seq
done

#Post process pile of sams.
echo "Merging and removing sams..."
samtools merge -f ${out}.sam ${readsDir}/*.sam
#rm $readsDir/*.sam
echo "Make sure to check and remove sams..."

echo "Making bam output..."
samtools view -S -b ${out}.sam > ${out}.bam
samtools sort -o ${out}_sort.bam ${out}.bam
samtools index ${out}_sort.bam
echo "Done. "
