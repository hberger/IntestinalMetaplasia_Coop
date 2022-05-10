BAM_FOLDER=../DeDup

CONTRA_FOLDER=/data_genome1/SharedSoftware/CONTRA/CONTRA.v2.0.8
#REGIONS=/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/S07604624_Covered_fixed_chroms.bed
#REF_FASTA=$(readlink -f /data_genome1/SharedSoftware/GATK/resources/human_g1k_v37_decoy.fasta)
REGIONS=/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed
REF_FASTA=/data_genome1/References/Human/Sequences/Genome/GRCh38_no_alt_plus_hs38d1_Verily/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

R_PATH=/data_genome1/SharedSoftware/R/4.1.0/bin
export R_LIBS_USER=~/R/packages/4.1

all_chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

control="4447_N" # 592
samples="4447_C 4447_D 4447_E 4447_F 4447_G 4447_H"

for sample in $samples; do
        ANALYSIS_FOLDER=./CN_${sample}

        PATH=$R_PATH:$PATH python ${CONTRA_FOLDER}/contra.py -t $REGIONS -s $BAM_FOLDER/${sample}_merged_dedup_dedup.bam -c $BAM_FOLDER/${control}_merged_dedup_dedup.bam -o $ANALYSIS_FOLDER -l

done

control="4447_O" # 594
samples="4447_I 4447_J 4447_K 4447_L 4447_P 4447_Q"

for sample in $samples; do
        ANALYSIS_FOLDER=./CN_${sample}

        PATH=$R_PATH:$PATH python ${CONTRA_FOLDER}/contra.py -t $REGIONS -s $BAM_FOLDER/${sample}_merged_dedup_dedup.bam -c $BAM_FOLDER/${control}_merged_dedup_dedup.bam -o $ANALYSIS_FOLDER -l

done

control="4447_M" # 591
samples="4447_B"


for sample in $samples; do
        ANALYSIS_FOLDER=./CN_${sample}

        PATH=$R_PATH:$PATH python ${CONTRA_FOLDER}/contra.py -t $REGIONS -s $BAM_FOLDER/${sample}_merged_dedup_dedup.bam -c $BAM_FOLDER/${control}_merged_dedup_dedup.bam -o $ANALYSIS_FOLDER -l

done
