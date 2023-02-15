import vcf
import pandas as pd
import argparse
import os
import re


def merge_vep_annovar(work_path, vep_path, annovar_path):
    '''
    Merge VEP and ANNOVAR annotation files
    :param work_path:
    :param vep_path:
    :param annovar_path:
    :return:
    '''

    os.system('mkdir {}/tmp'.format(work_path))
    bcftools = '/share/inspurStorage/home1/caojx/anaconda3/envs/detect_sv/bin/bcftools'

    # Loop the 24 chromosome
    for chr_ in list(range(1,23)) + ['X', 'Y']:
        # VEP annotation
        chr_vep_annot = os.path.join(vep_path, 'cohort_SNV_HC_chr{}_annot.vcf.gz'.format(chr_))

        # VEP annot fields
        os.system('{} +split-vep -l {} > {}/vep_annot_header.txt'.format(bcftools, chr_vep_annot, work_path))
        df_header = pd.read_csv('{}/vep_annot_header.txt'.format(work_path), sep='\t', names=['idx', 'annot'])
        dict_header = dict(zip(df_header['idx'], df_header['annot']))

        # select only first 8 columnes
        awk_args = '\"{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}\"'
        os.system("{} view -H {} | awk {} > {}/tmp/cohort_SNV_HC_chr{}_vep.txt".format(bcftools,chr_vep_annot,awk_args, work_path, chr_))

        # Process vep annotation files
        df_vep_1 = pd.read_csv('{}/tmp/cohort_SNV_HC_chr{}_vep.txt'.format(work_path, chr_),sep=' ',low_memory=False,header=None,names=['vcf_chr','vcf_pos','vcf_id','vcf_ref','vcf_alt','vcf_qual','vcf_filter','vcf_info'])
        df_vep_1['vcf_identifier'] = df_vep_1.apply(lambda x: str(x['vcf_chr']) + '-' + str(x['vcf_pos']) + '-' + str(x['vcf_ref']) + '-' + str(x['vcf_alt']),axis=1)
        df_vep_1['vep_impact'] = df_vep_1['vcf_info'].apply(lambda x: x.split(';')[-1].split('|')[2] if 'CSQ' in x else 'UNK')
        df_vep_1['vep_consequence'] = df_vep_1['vcf_info'].apply(lambda x: x.split(';')[-1].split('|')[1] if 'CSQ' in x else 'UNK')
        df_vep_1['vep_gene'] = df_vep_1['vcf_info'].apply(lambda x: x.split(';')[-1].split('|')[4] if 'CSQ' in x else 'UNK')
        df_vep_1['cohort_AF'] = df_vep_1['vcf_info'].apply(lambda x:re.findall('(?<=AF=).+?(?=;)',x)[0])

        # Annovar annotation
        chr_annovar_annot = os.path.join(annovar_path, 'cohort_SNV_HC_chr{}.hg19_multianno.txt'.format(chr_))

        # Process annovar annotation files
        df_annovar = pd.read_csv(chr_annovar_annot, sep='\t', low_memory=False)

        # rename the duplicates columns
        col_annovar_1 = df_annovar.columns.tolist()
        col_annovar_2 = [x.replace('.1','_genome') if ('.1' in x and 'AF' in x) else x for x in col_annovar_1]
        df_annovar.columns = col_annovar_2

        # Merges
        df_merge = pd.concat([df_vep_1,df_annovar],axis=1)
        df_merge.to_csv('{}/cohort_SNV_HC_chr{}_annotation.csv'.format(work_path, chr_),index=False)

        # Clean
        os.system('rm -rf {}/tmp'.format(work_path))



def filter_coding_variants(df_annot):
    '''
    Filter rare, likely deleterious coding variants
    :param df_annot:
    :return:
    '''

    # coding
    df_coding = df_annot[df_annot['vep_impact'] != 'MODIFIER'].copy()

    # rare
    df_coding_1 = df_coding[(df_coding['AF_eas'] < 0.01) & (df_coding['cohort_AF'] < 0.01)].copy()

    # likely pathogenic
    df_coding_1_high = df_coding_1[df_coding_1['vep_impact'] == 'HIGH'].copy()
    df_coding_1_high['patho_num'] = 10
    df_coding_1_middle = df_coding_1[df_coding_1['vep_impact'] != 'HIGH'].copy()
    df_coding_1_middle['patho_num'] = df_coding_1_middle.apply(lambda x: [x['SIFT_pred'] in ['P', 'D'], x['Polyphen2_HVAR_pred'] in ['P', 'D'],
                                                                          x['LRT_pred'] in ['D'],x['MutationTaster_pred'] in ['A', 'D'], x['MetaSVM_pred'] in ['D'],
                                                                          x['FATHMM_pred'] in ['D'],x['REVEL_score'] > 0.5, x['CADD_phred'] > 20].count(True), axis=1)
    df_coding_2 = pd.concat([df_coding_1_high, df_coding_1_middle])
    df_coding_3 = df_coding_2[(df_coding_2['patho_num'] >= 4) | (df_coding_2['InterVar_automated'].isin(['Likely pathogenic', 'Pathogenic']))].copy()

    return df_coding_3



def filter_noncoding_variants(df_annot):
    '''
    Filter rare, likely deleterious non-coding variants
    :param df_annot:
    :return:
    '''

    # noncoding
    df_noncoding = df_annot[df_annot['vep_impact'] == 'MODIFIER'].copy()

    # rare
    df_noncoding_1 = df_noncoding[(df_noncoding['AF_eas_genome'] < 0.01) & (df_noncoding['cohort_AF'] < 0.01)].copy()

    # likely pathogenic
    df_noncoding_1['patho_num'] = df_noncoding_1.apply(lambda x: [x['GWAVA_region_score'] > 0.5 or x['GWAVA_tss_score'] > 0.5, x['CADD13_PHRED'] > 15,
                                                                       x['FATHMM_noncoding'] > 0.9, x['dann'] > 0.9, x['Eigen'] > 1, x['gerp++gt2'] > 2].count(True),axis=1)
    df_noncoding_2 = df_noncoding_1[(df_noncoding_1['patho_num'] >= 3)].copy()

    return df_noncoding_2



def filter_rare_deleterious_variants(work_path, merge_annot_path):

    os.system('mkdir -p {}/table_format'.format(work_path))

    # Loop the 24 chromosome
    for chr_ in list(range(1, 23)) + ['X', 'Y']:

        # Annot files
        df_annot_1 = pd.read_csv('{}/cohort_SNV_HC_chr{}_annotation.csv'.format(merge_annot_path, chr_),low_memory=False)

        # Change columns dtype
        col_AF_exome = ['AF','AF_popmax','AF_male','AF_female','AF_raw','AF_afr','AF_sas','AF_amr','AF_eas','AF_nfe','AF_fin','AF_asj','AF_oth']
        col_AF_genome = ['AF_genome','AF_popmax_genome','AF_male_genome','AF_female_genome','AF_raw_genome','AF_afr_genome','AF_sas_genome','AF_amr_genome','AF_eas_genome','AF_nfe_genome','AF_fin_genome','AF_asj_genome','AF_oth_genome']
        col_patho_genome = ['dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','dpsi_max_tissue','dpsi_zscore','GWAVA_region_score','GWAVA_tss_score','GWAVA_unmatched_score','gerp++gt2','CADD13_PHRED','FATHMM_noncoding','dann','Eigen']
        col_patho_exome = ['CADD_phred','REVEL_score']
        df_annot_1.loc[:,col_AF_exome] = df_annot_1.loc[:,col_AF_exome].applymap(lambda x:0 if x == '.' else float(x))
        df_annot_1.loc[:,col_AF_genome] = df_annot_1.loc[:,col_AF_genome].applymap(lambda x:0 if x == '.' else float(x))
        df_annot_1.loc[:,col_patho_genome] = df_annot_1.loc[:,col_patho_genome].applymap(lambda x:-100 if x == '.' else float(x))
        df_annot_1.loc[:,col_patho_exome] = df_annot_1.loc[:,col_patho_exome].applymap(lambda x:-100 if x == '.' else float(x))
        df_annot_1['cohort_AF'] = df_annot_1['cohort_AF'].fillna(0).astype(float)

        # Filter coding variants
        df_coding_filtered = filter_coding_variants(df_annot_1)
        df_coding_filtered.to_csv('{}/table_format/cohort_SNV_HC_chr{}_coding_rp.csv'.format(work_path,chr_), index=False)

        # Filter noncoding variants
        df_noncoding_filtered = filter_noncoding_variants(df_annot_1)
        df_noncoding_filtered.to_csv('{}/table_format/cohort_SNV_HC_chr{}_noncoding_rp.csv'.format(work_path,chr_), index=False)

    # Merged
    df_filtered_annotation = pd.DataFrame()
    for line_ in os.listdir('{}/table_format'.format(work_path)):
        df_line_ = pd.read_csv('{}/table_format/{}'.format(work_path,line_), low_memory=False)
        df_filtered_annotation = pd.concat([df_filtered_annotation, df_line_])
    df_filtered_annotation.to_csv("{}/cohort_SNV_HC_rp.csv".format(work_path), index=False)

    return df_filtered_annotation




def generate_vcf(work_path, vep_path):

    os.system('mkdir -p {}/vcf_format'.format(work_path))
    os.system('mkdir -p {}/vcf_format/tmp'.format(work_path))

    bcftools = '/share/inspurStorage/home1/caojx/anaconda3/envs/detect_sv/bin/bcftools'

    # vcf_header
    os.system("{} view -h {}/cohort_SNV_HC_chr1_annot.vcf.gz | grep -v '#CHROM' > {}/vcf_format/tmp/vcf_header.txt".format(bcftools, vep_path, work_path))

    # Loop the 24 chromosome
    for chr_ in list(range(1, 23)) + ['X', 'Y']:
        df_coding_filtered = pd.read_csv('{}/table_format/cohort_SNV_HC_chr{}_coding_rp.csv'.format(work_path,chr_), low_memory=False)
        df_noncoding_filtered = pd.read_csv('{}/table_format/cohort_SNV_HC_chr{}_noncoding_rp.csv'.format(work_path,chr_), low_memory=False)
        df_filtered = pd.concat([df_coding_filtered, df_noncoding_filtered])
        df_filtered_1 = df_filtered.iloc[:,:7]
        df_filtered_1.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER']
        df_filtered_1.to_csv('{}/vcf_format/tmp/{}_filtered.txt'.format(work_path, chr_),sep='\t',index=False)

        # make a vcf
        os.system("cat {work_path}/vcf_format/tmp/vcf_header.txt {work_path}/vcf_format/tmp/{chr}_filtered.txt > {work_path}/vcf_format/tmp/{chr}_filtered.vcf".format(work_path=work_path, chr=chr_))
        os.system("{bcftools} sort {work_path}/vcf_format/tmp/{chr}_filtered.vcf -Oz -o {work_path}/vcf_format/tmp/{chr}_filtered.vcf.gz".format(bcftools=bcftools, work_path=work_path, chr=chr_))
        os.system("tabix {work_path}/vcf_format/tmp/{chr}_filtered.vcf.gz".format(work_path=work_path, chr=chr_))

        # extract those rare and likely deleterious variants from raw vcf
        os.system("bcftools isec -p {work_path}/vcf_format/{chr} -n=2 {work_path}/vcf_format/tmp/{chr}_filtered.vcf.gz {vep_path}/cohort_SNV_HC_chr{chr}_annot.vcf.gz".format(work_path=work_path, vep_path=vep_path, chr=chr_))
        os.system("{bcftools} sort {work_path}/vcf_format/{chr}/0001.vcf -Oz -o {work_path}/vcf_format/cohort_SNV_HC_chr{chr}_rp.vcf.gz".format(bcftools=bcftools, work_path=work_path, chr=chr_))
        os.system("tabix {work_path}/vcf_format/cohort_SNV_HC_chr{chr}_rp.vcf.gz".format(work_path=work_path, chr=chr_))

    os.system("find {work_path}/vcf_format | grep _rp.vcf.gz  | grep -v .tbi > {work_path}/vcf_format/tmp/vcf_list.txt".format(work_path=work_path))
    os.system("bcftools concat -f {work_path}/vcf_format/tmp/vcf_list.txt -a -Oz -o {work_path}/vcf_format/merged_rare_pathogenic.vcf.gz".format(work_path=work_path))
    os.system("bcftools sort {work_path}/vcf_format/merged_rare_pathogenic.vcf.gz -Oz -o {work_path}/merged_sorted_rare_pathogenic.vcf.gz".format(work_path=work_path))
    os.system("tabix {work_path}/merged_sorted_rare_pathogenic.vcf.gz".format(work_path=work_path))

