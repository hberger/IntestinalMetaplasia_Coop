if [ $# -lt 3 ]; then
   echo "Usage:"
   echo "DNASeq_calling_pipeline_single_sample BAM_FOLDER SAMPLE ROI.bed [MAX_COVERAGE]"
   exit
fi


BAM_FOLDER=$1
sample=$2
ROI=$3

if [ $# -eq 4 ]; then
   MAX_COV=$4
else
   MAX_COV=6000
fi

echo "MAX_COV: " $MAX_COV

GATK=/data_genome1/SharedSoftware/GATK/gatk-4.1.9.0/gatk
PICARD_FOLDER=/data_genome1/SharedSoftware/Picard
REF_FASTA=../../../mapping/DNA_h38/GenomeIndices/GRCh38_no_alt_plus_hs38d1_Verily_sequence/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
KNOWN_SITES=/data_genome1/SharedSoftware/GATK/resources_hg38/dbsnp_146.hg38.vcf

# put the temp folder where we know there will be enough space available 
TMP_FOLDER=$(pwd)/tmp
if ! [ -e $TMP_FOLDER ]; then mkdir $TMP_FOLDER; fi

#BASH specific
#http://stackoverflow.com/questions/1527049/bash-join-elements-of-an-array
function join { local IFS="$1"; shift; echo "$*"; }

	all_sample_lanes=""
        all_sample_lane_cnt=0

    for sample_lane in ${BAM_FOLDER}/${sample}*_sorted.bam; do

	INPUT_BAM_FILE=$sample_lane
	i=$(basename $sample_lane _sorted.bam)

	run_lane=$(echo $i | awk 'BEGIN{FS="_"} {print $3"_"$4;}')
	RG=${run_lane}_${sample}
	LIBRARY=$sample
	PLATFORM="illumina"
	RGPU=$run_lane
	SAMPLE=$sample

	RG_ADDED_SORTED_BAM=${i}_rg.bam
	if ! [ -e $RG_ADDED_SORTED_BAM ]; then
	   samtools view -F 8 -u -O BAM $INPUT_BAM_FILE | samtools addreplacerg -r "@RG\tID:${RG}\tPL:${PLATFORM}\tPU:${RGPU}\tLB:${LIBRARY}\tSM:${SAMPLE}" -@ 4 -O BAM -o $RG_ADDED_SORTED_BAM - 
	   #java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/AddOrReplaceReadGroups.jar I=$INPUT_BAM_FILE O=$RG_ADDED_SORTED_BAM SO=coordinate RGID=$RG RGLB=$LIBRARY RGPL=$PLATFORM RGPU=$RGPU RGSM=$sample TMP_DIR=$TMP_FOLDER
	fi

	DEDUPPED_BAM=${i}_rg_dedup.bam
	DEDUP_METRICS_FILE=${i}.dedup.metrics
        if ! [ -e $DEDUPPED_BAM ]; then
	    java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=$RG_ADDED_SORTED_BAM O=$DEDUPPED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER
	    #$GATK MarkDuplicatesSpark -I $RG_ADDED_SORTED_BAM -O $DEDUPPED_BAM --conf 'spark.executor.cores=1' --tmp-dir $TMP_FOLDER
            #samtools -@ 4 index $DEDUPPED_BAM
        fi


        RECAL_TAB_FILE=${i}_recal_data.tab
        POST_RECAL_TAB_FILE=${i}_post_recal_data.tab
        RECAL_PLOTS=${i}_recalibration_plots.pdf
        RECALIBRATED_BAM=${i}_rg_dedup_recal.bam
        
        SPARK_WORK_PATH=./tmp
        if ! [ -e $RECALIBRATED_BAM ]; then
        
           $GATK --java-options -Xmx8G \
             BaseRecalibratorSpark \
             -R $REF_FASTA -I $DEDUPPED_BAM --known-sites $KNOWN_SITES \
             -O $RECAL_TAB_FILE -L $ROI \
             -- --spark-runner LOCAL --spark-master local[12] \
             --conf spark.local.dir=$SPARK_WORK_PATH

           $GATK --java-options -Xmx8G \
             ApplyBQSRSpark -R $REF_FASTA \
             -I $DEDUPPED_BAM -bqsr $RECAL_TAB_FILE -L $ROI \
             --static-quantized-quals 10 --static-quantized-quals 20 \
             --static-quantized-quals 30 -O $RECALIBRATED_BAM \
             -- --spark-runner LOCAL --spark-master local[12] \
             --conf spark.local.dir=$SPARK_WORK_PATH
        
	   #$GATK BaseRecalibrator -I $DEDUPPED_BAM -R $REF_FASTA --known-sites $KNOWN_SITES -O $RECAL_TAB_FILE
           #$GATK ApplyBQSR -R $REF_FASTA -I $DEDUPPED_BAM --bqsr-recal-file $RECAL_TAB_FILE -O $RECALIBRATED_BAM
        fi
        
        if [ -e $RECALIBRATED_BAM ]; then
            #only remove intermediate files if we reached the last step
            rm $RG_ADDED_SORTED_BAM #$DEDUPPED_BAM
        fi
#	rm $RG_ADDED_SORTED_BAM        

	all_sample_lanes+="INPUT="
	all_sample_lanes+=$i
	all_sample_lanes+="_rg_dedup_recal.bam "
        all_sample_lane_cnt+=1

   done 

   # now merge per-lane BAMs to per-sample BAMs
   #inputfiles=$(join ' INPUT=' $all_sample_lanes)
   echo $all_sample_lanes
   inputfiles=$all_sample_lanes

   MERGED_BAM=${sample}_merged.bam
   if [ $all_sample_lane_cnt -gt 1 ]; then 

	   if ! [ -e $MERGED_BAM ]; then
	       java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MergeSamFiles.jar $inputfiles OUTPUT=$MERGED_BAM TMP_DIR=$TMP_FOLDER 
	   fi

	   # optional: rerun dedup and realign on per-sample BAMs
	   #DEDUPPED_MERGED_BAM=${sample}_merged_dedup.bam
	   #DEDUP_METRICS_FILE=${sample}.merged.dedup.metrics
	   #if ! [ -e $DEDUPPED_MERGED_BAM ]; then
	   #    java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=${sample}_merged.bam O=$DEDUPPED_MERGED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER 
	   #fi


   else
	REALIGNED_MERGED_BAM=${sample}_merged_dedup.bam
	REALIGNED_MERGED_BAI=${sample}_merged_dedup.bai
	cp $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g') $REALIGNED_MERGED_BAM
	cp $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g' | sed -e 's/.bam/.bam.bai/') $REALIGNED_MERGED_BAI
   fi

 #  if [ -e $DEDUPPED_MERGED_BAM ]; then
 #      rm ${sample}*_rg_dedup.bam ${sample}*_rg.bam
 #  fi

