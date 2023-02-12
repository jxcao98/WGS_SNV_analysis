#!/bin/bash

#################################################################################

## Step4: Annotation
## Usage: ./run_04_Annotation.sh myConfigFile.config

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
runDir=${outputFolderPath}/02_SNV/03_Annotation


# =============================================================================
#                           Run Sentieon Pipeline
# =============================================================================
# [0] prepare input files
## Split vcf by chromosomes
# We recommand split a cohort vcf by chromosomes to accelerate annotation.
mkdir -p ${runDir}/00_input
mkdir -p ${runDir}/00_input/vcf_by_chr
cd ${runDir}/00_input/vcf_by_chr

# get chromosomes names
tabix --list-chroms ${outputFolderPath}/02_SNV/01_joint_calling/02_vqsr/cohort_SNV_HC.vcf.gz > chromosomes_all.txt
head -24 chromosomes_all.txt > chromosomes.txt

# split vcf
while IFS= read -r line; do
  tabix -h ${outputFolderPath}/02_SNV/01_joint_calling/02_vqsr/cohort_SNV_HC.vcf.gz $line | bgzip -c > cohort_SNV_HC_chr$line.vcf.gz
  tabix cohort_SNV_HC_chr$line.vcf.gz;
done < chromosomes.txt


## Change vcf to annovar.input format
mkdir -p ${runDir}/00_input/vcf2annovarinput
cd ${runDir}/00_input/vcf2annovarinput

for vcf_chr in $(ls ${runDir}/00_input/vcf_by_chr)
do
	array=(${vcf_chr//./ })
	new_file_name=${array[0]}
	$ANNOVAR/convert2annovar.pl -format vcf4 ${runDir}/00_input/vcf_by_chr/$vcf_chr -outfile $new_file_name.avinput -allsample -includeinfo -withfreq
	# Select only the first 5 columns to run faster
	cat $new_file_name.avinput | awk '{print $1,$2,$3,$4,$5}' > ${new_file_name}_5cols.avinput
done



# [1] VEP
mkdir ${runDir}/vep
mkdir ${runDir}/vep/log

cd ${runDir}/vep

# Loop over each chromosome
vcf_path=${runDir}/00_input/vcf_by_chr
for vcf_chr in $(ls $vcf_path | grep vcf.gz | grep -v tbi)
do	
{

	array=(${vcf_chr//./ })
	file_name=${array[0]}
	
	logfile=./log/$file_name.log
	exec >$logfile 2>&1

	vep --cache --offline --format vcf --vcf --assembly GRCh37 --force_overwrite \
		--dir_cache $VEP_PATH \
		--dir_plugins $VEP_PATH/Plugins \
		--input_file $vcf_path/$file \
		-o stdout \
		--fasta $VEP_PATH/homo_sapiens/108_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
		--sf ./log/$file_name.summary.html \
		--warning_file ./log/$file_name.warning.txt \
		--everything \
		--fork $nt \
		--pick \
		--pick_order canonical,biotype,rank,length \
		--plugin LoF,loftee_path:$VEP_PATH/Plugins,human_ancestor_fa:$VEP_PATH/Plugins/human_ancestor_fa/human_ancestor.fa.gz,phylocsf_data:$VEP_PATH/Plugins/phylocsf_gerp.sql.gz \
		--plugin dbNSFP,$VEP_PATH/dbnsfp4/dbNSFP4.1a_grch37.gz,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,CADD_phred,FATHMM_pred,MetaSVM_pred,REVEL_score \
		| bgzip -c > ${file_name}_annot.vcf.gz
	tabix ${file_name}_annot.vcf.gz

}
done




# [2] ANNOVAR
mkdir ${runDir}/annovar
mkdir ${runDir}/annovar/log

cd ${runDir}/annovar

# Loop over each chromosome
annovarinput_path=${runDir}/00_input/vcf2annovarinput
for file_chr in $(ls $annovarinput_path | grep _5cols)
do	
{

	array=(${file_chr//_5cols/ })
	file_name=${array[0]}
	
	logfile=./log/$file_name.log
	exec >$logfile 2>&1

	$ANNOVAR_PATH/table_annovar.pl \
		$annovarinput_path/$file_chr \
		$ANNOVAR_PATH/humandb/ \
		--buildver hg19 \
		--out ./$file_name \
		--remove \
		--protocol ensGene,cytoBand,intervar_20180118,dbnsfp41a,dbscsnv11,spidex,regsnpintron,gwava,gerp++gt2,cadd13,fathmm,eigen,dann,clinvar_20210501,gnomad211_exome,gnomad211_genome,gwasCatalog,tfbsConsSites \
		--operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r \
		--nastring . \
		--thread $nt

}
done
