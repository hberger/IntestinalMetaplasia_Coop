#!/usr/bin/bash
# Prerequisites: 
# - download human and mouse genome sequences and prepare the genome index (see subfolders in GenomeIndices)
# - Have python and package pysam installed (pip install pysam)

bash ./BWA_MEM_with_mm38/align_BWA_MEM.sh
bash ./BWA_MEM_split/SCRIPT_split
bash ./BWA_MEM_hs_only/SCRIPT
bash ./BWA_MEM_hs_only_merged/merge_BAMs.sh
