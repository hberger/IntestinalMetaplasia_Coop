BAM_FOLDER=../GATK/processed_BAMs

samples=$(ls ${BAM_FOLDER}/*.bam)

for s in $samples; do 
	samtools rmdup -s $s $(basename $s .bam)_dedup.bam
	samtools index $(basename $s .bam)_dedup.bam
done

