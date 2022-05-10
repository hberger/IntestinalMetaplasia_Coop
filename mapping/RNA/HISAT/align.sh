#!/bin/bash
#ulimit -v 30258917744
HISAT_FOLDER=/data_genome1/SharedSoftware/HISAT/download/hisat2-2.1.0
#GENOME_DIR=/data_genome1/References/GenomeIndices/Human/HISAT2/hg38/genome 
GENOME_DIR=/data_genome1/References/GenomeIndices/Human/HISAT2/GRCh38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set
SPLICE_SITES=/data_genome1/References/Human/Annotations/hg38/Ensembl/splice_sites.txt

FASTQ_PATH=../../../seqs/FASTQ

export PATH=$PATH:$HISAT_FOLDER

function join { local IFS="$1"; shift; echo "$*"; }

TMP_FOLDER=./tmp

mkdir $TMP_FOLDER

tail -n +2 ../../../metadata/RNA/SraRunTable.txt | cut -d "," -f 1 | while read sample_line
do # 1st pass
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=${FASTQ_PATH}/${sample}_1.fastq.gz
  reads2=${FASTQ_PATH}/${sample}_2.fastq.gz
#  $STARexec --genomeDir $GENOME_DIR --readFilesIn $reads1 --readFilesCommand zcat --outFileNamePrefix ${sample}_ --runThreadN 8 --outSAMstrandField intronMotif --genomeLoad NoSharedMemory --twopassMode Basic

  $HISAT_FOLDER/hisat2 -p 12 --known-splicesite-infile $SPLICE_SITES  -x $GENOME_DIR -1 $reads1 -2 $reads2 -S ${sample}_Aligned.out.sam

  /data_genome1/SharedSoftware/samtools_htslib/samtools view -u -b ${sample}_Aligned.out.sam | /data_genome1/SharedSoftware/samtools_htslib/samtools sort -@ 4 -T ${TMP_FOLDER}/CHUNK -O bam -o ${sample}_sorted.bam
  /data_genome1/SharedSoftware/samtools_htslib/samtools index -b ${sample}_sorted.bam

  rm ${sample}_Aligned.out.sam

done

