#!/bin/bash

pipeline=./

# [Begin] Make a config file, and a commas seperated file containing cohort infomation (sample_name, fq1_path, fq2_path)
# Example: ./cohort.config, ./cohort_info.txt

# [0] Fastqc
sh $pipeline/Step0_run_fastqc.sh cohort.config

# [1] Run sentieon on per sample (Sample Level)
sh $pipeline/Step1_run_sentieon.sh cohort.config

# [2] Joint calling and VQSR (Cohort Level)
sh $pipeline/Step2_joint_calling.sh cohort.config
# You need to change the plink.fam files manually according to the cohort information.

# [3] QC
sh $pipeline/Step3_QC.sh cohort.config

# [4] Annotation
sh Step4_Annotation.sh cohort.config
