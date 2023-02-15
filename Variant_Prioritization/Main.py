"""
Filtering to keep rare and likely deleterious variants according to annotation results.
A typical VCF file and a view-friendly table file containing the annotation information are generated.
"""

import pandas as pd
import os
import numpy as np
import Variant_filter
import Parse_vcf


# [1] Merge Annovar and VEP annotation files
vep_annot_path = 'outputFolderPath/02_SNV/03_Annotation/vep'
annovar_annot_path = 'outputFolderPath/02_SNV/03_Annotation/annovar'
Step1_work_path = 'outputFolderPath/02_SNV/03_Annotation/merge'

Variant_filter.merge_vep_annovar(Step1_work_path, vep_annot_path, annovar_annot_path)


# [2] Filter variants (keep rare, likely deleterious variants)
Step2_work_path = 'outputFolderPath/02_SNV/04_rare_variants'
Step2_run_path = 'outputFolderPath/02_SNV/04_rare_variants/01_filter_variants'
os.system('mkdir -p {}'.format(Step2_work_path))
os.system('mkdir -p {}'.format(Step2_run_path))

df_filter_annotation = Variant_filter.filter_rare_deleterious_variants(Step2_run_path, Step1_work_path)
Variant_filter.generate_vcf(Step2_run_path, vep_annot_path)


# [3] Parse Vcf to get cohort information
Step3_run_path = 'outputFolderPath/02_SNV/04_rare_variants/02_parse_vcf'
ped_file = 'outputFolderPath/02_SNV/01_joint_calling/03_vcf2plink/cohort_SNV_HC_v3.fam'
samplesToRemove = 'outputFolderPath/02_SNV/03_QC/samplesToRemove.txt'
df_filter_cohort_info = Parse_vcf.vcf2table(os.path.join(Step2_run_path,'vcf_format/merged_sorted_rare_pathogenic.vcf.gz'), Step3_run_path, ped_file, samplesToRemove)


# [4] Merge
Step4_run_path = 'outputFolderPath/02_SNV/04_rare_variants/03_merge'
df_merge = pd.merge(df_filter_annotation, df_filter_cohort_info, on='vcf_identifier')


# By now, we have a set that contains only rare and likely deleterious variants.
# The next analyses are personalized, and examine loci in known known disease causative genes may be a good start.

# [5] Examining variants in known disease causative genes
known_genes = []
df_filter_0 = df_merge[df_merge['SYMBOL'].isin(known_genes)].copy() # variants in known disease causative genes
df_filter_1 = df_filter_0[(df_filter_0['NC_hom_num'] == 0) & (df_filter_0['NC_het_num'] == 0)].copy() # variants not carried by NC
df_filter_2 = df_filter_1[(df_filter_1['Case_hom_num'] > 0) | (df_filter_1['Case_het_num'] > 0)].copy() # Variants carried by Case

# We can also explore novel risk genes, and those loci that recurrent in cases may be promising.