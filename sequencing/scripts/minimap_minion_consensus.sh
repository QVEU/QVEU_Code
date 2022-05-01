#!/usr/bin/env bash
#minimap_nanopore_consensus.sh
#usage: sh minimap_nanopore_consensus.sh templatefasta.fa path-to-"fastq_pass"-Directory/
  #ARG. 1: template fasta file
  #ARG. 2: directory holding fastq files (usually this is "fastq_pass" directory in the minion run directory)

echo ${1}
echo ${2}

#Map reads in fastq files  and generate sam file.
echo "Mapping reads..."
for i in ${2}/*.fastq.gz
do
  echo $i
  /Users/dolanpt/opt/anaconda3/pkgs/minimap2-2.24-h1f540d2_1/bin/minimap2 -ax map-ont $1 $i > ${i/\.fastq.gz/}.sam      # for Oxford Nanopore reads
done

#Post process pile of sams.
echo "Merging and removing sams..."
cat ${2}/*.sam > ${2}/merge.sam
#samtools merge -f ${2}_out.sam ${2}/FAT*.sam
#rm ${2}/FAT*.sam

#picard MergeSamFiles -I  FAT*.sam -O  merged_files.bam

#make bam files
echo "Making bam output..."

samtools view -S -b ${2}/merge.sam > ${2}/merge.bam

#sort and pile up reads. make bed cov file with bed files from template.

samtools sort ${2}/merge.bam > ${2}/merge_sort.bam
samtools mpileup ${2}/merge_sort.bam > ${2}/merge_sort_pile.pile
#samtools bedcov -Q 8 MT276325.1.bed ${2}.sam
echo "Done. "
