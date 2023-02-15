import vcf
import os
import pandas as pd


def vcf2table(vcf_path, work_path, ped_file, samplesToRemove):

    os.system('mkdir -p {}/tmp'.format(work_path))

    # read vcf file
    vcf_reader = vcf.Reader(filename=vcf_path)
    file_name = vcf_path.split('/')[-1].split('.vcf')[0]

    # read ped file
    ped_file = pd.read_csv(ped_file, sep='\t', names = ['FID', 'IID', 'FAT', 'MOT', 'Gender', 'Diag'])
    ped_file['PHE'] = ped_file['Diag'].map({0:'Unknown', 1:'NC', 2:'Case'})

    # blacklist
    exclude_sample = pd.read_csv(samplesToRemove,sep='\t',names=['IID'])
    exclude_sample = exclude_sample['IID'].tolist()

    # get vep fields
    os.system('bcftools +split-vep -l {} > {}.header.txt'.format(vcf_path, work_path))
    df_header = pd.read_csv('{}.header.txt'.format(work_path), sep='\t', names=['idx', 'annot'])
    header_list = df_header['annot'].tolist()

    # generate information
    df_vcf = pd.DataFrame()
    n_row = 0

    # parse vcf
    for record in vcf_reader:
        n_row += 1
        if n_row % 1000 == 0:
            df_vcf.to_csv('{}/tmp/{}_line_{}.csv'.format(work_path, file_name, str(n_row)), index=False)
            df_vcf = pd.DataFrame()

        # base info
        df_vcf.loc[n_row, 'var_chrom'] = record.CHROM
        df_vcf.loc[n_row, 'var_pos'] = int(record.POS)
        df_vcf.loc[n_row, 'var_id'] = record.ID
        df_vcf.loc[n_row, 'var_ref'] = record.REF
        df_vcf.loc[n_row, 'var_alt'] = record.ALT[0]
        df_vcf.loc[n_row, 'var_qual'] = record.QUAL
        df_vcf.loc[n_row, 'var_num_called'] = record.num_called
        df_vcf.loc[n_row, 'var_call_rate'] = record.call_rate
        df_vcf.loc[n_row, 'var_num_unknown'] = record.num_unknown
        df_vcf.loc[n_row, 'var_num_hom_ref'] = record.num_hom_ref
        df_vcf.loc[n_row, 'var_num_het'] = record.num_het
        df_vcf.loc[n_row, 'var_num_hom_alt'] = record.num_hom_alt
        df_vcf.loc[n_row, 'vcf_identifier'] = str(record.CHROM) + '-' + str(int(record.POS)) + '-' + str(record.REF) + '-' + str(record.ALT[0])

        try:
            df_vcf.loc[n_row, header_list] = record.INFO['CSQ'][0].split('|')  # vep annotation fields
        except:
            df_vcf.loc[n_row, header_list] = -1

        het_sam_info = record.get_hets()
        hom_sam_info = record.get_hom_alts()
        miss_sam_info = record.get_unknowns()
        dict_sam_info = {'het': het_sam_info, 'hom': hom_sam_info, 'miss': miss_sam_info}

        for gt_ in dict_sam_info.keys():
            dict_phe = {}
            for phe_ in ped_file['PHE'].value_counts().index.tolist():
                dict_phe[phe_] = []

            # classify samples by phenotype
            for sam_info_ in dict_sam_info[gt_]:
                sam_ = sam_info_.sample
                if sam_ not in exclude_sample:
                    sam_phe_ = ped_file.loc[sam_, 'Diag']
                    GT = sam_info_['GT']
                    GQ = sam_info_['GQ']
                    dict_phe[sam_phe_].append(sam_ + ':' + str(GT) + ':' + str(GQ))  # name + GT + GQ

            # write all category
            for phe_ in ped_file['PHE'].value_counts().index.tolist():
                df_vcf.loc[n_row, phe_ + '_' + gt_] = ','.join(dict_phe[phe_])
                df_vcf.loc[n_row, phe_ + '_' + gt_ + '_num'] = len(dict_phe[phe_])

    df_vcf.to_csv('{}/tmp/{}_line_tail.csv'.format(work_path, file_name), index=False)

    # merged
    df_filtered_cohort_info = pd.DataFrame()
    for line_ in os.listdir('{}/tmp'.format(work_path)):
        df_line_ = pd.read_csv('{}/tmp/{}'.format(work_path,line_), low_memory=False)
        df_filtered_cohort_info = pd.concat([df_filtered_cohort_info, df_line_])
    df_filtered_cohort_info.to_csv("{}/rare_pathogenic_cohort_info.csv".format(work_path), index=False)

    return df_filtered_cohort_info