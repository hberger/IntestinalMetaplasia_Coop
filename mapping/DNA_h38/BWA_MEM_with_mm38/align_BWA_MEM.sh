#!/bin/bash

GENOME_INDEX=../GenomeIndices/GRCh38_noalt_hs38d_mm38/GRCh38_noalt_hs38d1_mm38.fasta

SAMTOOLS=/data_genome1/SharedSoftware/samtools_htslib/samtools
METAFILE_FOLDER=../../../metadata/DNA_both_batches
function join { local IFS="$1"; shift; echo "$*"; }


TMP_FOLDER=./samtools/BWA
mkdir -p $TMP_FOLDER
TMP_READS1_FIFO=${TMP_FOLDER}/tmp_reads1
TMP_READS2_FIFO=${TMP_FOLDER}/tmp_reads2

[ ! -p tmp_reads1 ] && mkfifo ${TMP_READS1_FIFO}
[ ! -p tmp_reads2 ] && mkfifo ${TMP_READS2_FIFO}

cat $METAFILE_FOLDER/files_by_sample_and_lane.txt | grep -v "IndexD710" | while read sample_line
do
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=$(echo $2 | sed -e 's/\,/ /g')
  reads2=$(echo $3 | sed -e 's/\,/ /g')


[ -e ${sample}_sorted.bam ] && continue

  zcat $reads1 > ${TMP_READS1_FIFO} &
  zcat $reads2 > ${TMP_READS2_FIFO} &

  bwa mem -t 24 $GENOME_INDEX ${TMP_READS1_FIFO} ${TMP_READS2_FIFO} | $SAMTOOLS view -u -b - | $SAMTOOLS sort -@ 4 -T ${TMP_FOLDER}/CHUNK -O bam -o ${sample}_sorted.bam
  $SAMTOOLS index -b ${sample}_sorted.bam

done

rm ${TMP_READS1_FIFO}
rm ${TMP_READS2_FIFO}

