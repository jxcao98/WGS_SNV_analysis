#!/bin/bash

#################################################################################

## Step4: Annotation
## Usage: ./Step4_Annotation.sh myConfigFile.config

#################################################################################

# =============================================================================
#                           Load Config Files
# =============================================================================

display_usage() {
        echo "Script must run with the configuration file."
        echo -e "\nUsage:\n $0 shortReadPipeline.config \n"
        }
# if less than one argument supplied, display usage
if [[ $# -eq 0 ]] ; then
    display_usage
    exit 1
fi
# check whether user had supplied -h or --help . If yes display usage
if [[ ( $# == "--help") ||  $# == "-h" ]] ; then
    display_usage
    exit 0
fi
# load configFile
source $1
if ! [ -n "$FastqList" ]; then
    echo "Something is wrong. FastqList var is unset. Please check if the config file is OK.";
    exit 1
else
    echo "Loading config file seems OK. FastqList='$FastqList'.";
    lastline=$(tail -n1 $1)
    if [ "$lastline" = "#### End of file ####" ]; then
        echo "Last line seems OK.";
    else
        echo "Last line must be '#### End of file ####'. Please check if the config file is OK.";
        exit 1
    fi
fi


# =============================================================================
#                           Set folders and parameters
# =============================================================================
mkdir -p ${outputFolderPath} || exit 1
mkdir -p ${outputFolderPath}/02_SNV || exit 1
mkdir -p ${outputFolderPath}/02_SNV/03_Annotation || exit 1
RunFolderPath=${outputFolderPath}/02_SNV/03_Annotation


# =============================================================================
#                           Run Sentieon Pipeline
# =============================================================================
mkdir ${RunFolderPath}/01_vep
mkdir ${RunFolderPath}/01_vep/01_vep_annot
mkdir ${RunFolderPath}/01_vep/01_vep_annot/log
vcf_path=${outputFolderPath}/02_SNV/01_joint_calling/03_change_format/01_vcf_by_chr


cd ${RunFolderPath}/01_vep/01_vep_annot

for file in $(ls $vcf_path | grep vcf.gz | grep -v tbi)
do	
{

	array=(${file//./ })
	new_file_name=${array[0]}
	
	logfile=${RunFolderPath}/01_vep/01_vep_annot/log/$new_file_name.log
	exec >$logfile 2>&1

	vep --cache --offline --format vcf --vcf --assembly GRCh37 --force_overwrite \
		--dir_cache $vep_path \
		--dir_plugins $vep_path/Plugins \
		--input_file $vcf_path/$file \
		-o stdout \
		--fasta $vep_path/homo_sapiens/108_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
		--sf ${RunFolderPath}/01_vep/01_vep_annot/log/$new_file_name.summary.html \
		--warning_file ${RunFolderPath}/01_vep/01_vep_annot/log/$new_file_name.warning.txt \
		--everything \
		--fork 20 \
		--pick \
		--pick_order canonical,biotype,rank,length \
		--plugin LoF,loftee_path:$vep_path/Plugins,human_ancestor_fa:$vep_path/Plugins/human_ancestor_fa/human_ancestor.fa.gz,phylocsf_data:$vep_path/Plugins/phylocsf_gerp.sql.gz \
		--plugin dbNSFP,$vep_path/dbnsfp4/dbNSFP4.1a_grch37.gz,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,CADD_phred,FATHMM_pred,MetaSVM_pred,REVEL_score \
		| bgzip -c > ${new_file_name}_annot.vcf.gz
	tabix ${new_file_name}_annot.vcf.gz

}
done
