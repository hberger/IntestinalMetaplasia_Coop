#!/bin/bash

# adapt that script to your needs. 
# Goal is to get a single BAM file + index for each sample
# If there are several lanes/FASTQ files per input sample, add those sample IDs to "samples" below and run samtools merge on them (lower part)
# If there are only single lanes/FASTQ files per input sample, you can use a symlink right away (upper part) 

SAMTOOLS=/data_genome1/SharedSoftware/samtools_htslib/samtools
METAFILE_FOLDER=../../../metadata/DNA_both_batches
function join { local IFS="$1"; shift; echo "$*"; }

samples="4447_M 4447_P"

cat $METAFILE_FOLDER/all_samples_aggregated.txt | while read sample_line
do
   [ -z "$sample_line" ] && continue
   set $sample_line
   sample=$(echo $1)

   ln -s ../BWA_MEM_hs_only/${sample}*.bam ${sample}_merged_sorted.bam
   ln -s ../BWA_MEM_hs_only/${sample}*.bai ${sample}_merged_sorted.bai

done

for sample in $samples
do
#  [ -z "$sample_line" ] && continue
#  set $sample_line
#  sample=$(echo $1)

  rm ${sample}_merged_sorted.*

  OUTFILE=${sample}_merged_sorted.bam

  $SAMTOOLS merge -@ 4 -O BAM $OUTFILE ../BWA_MEM_hs_only/${sample}_*.bam
  $SAMTOOLS index -b $OUTFILE

done

