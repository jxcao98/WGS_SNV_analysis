#!/bin/bash

#################################################################################

## Step2: Cohort-level SNV and QC
## Usage: ./Step2_joint_calling.sh myConfigFile.config

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

# Sentieon License
export SENTIEON_LICENSE=${SENTIEON_LICENSE}

# log
logfile=${outputFolderPath}/02_SNV/cohort_level_snv_detect.log
exec >$logfile 2>&1

# =============================================================================
#                           Run Pipeline
# =============================================================================

### [1] Joint Calling ###
mkdir -p ${outputFolderPath}/02_SNV/01_joint_calling || exit 1
cd ${outputFolderPath}/02_SNV/01_joint_calling

# Joint calling
gvcf_list=`ls ${outputFolderPath}/01_sentieon/*/*/06_call_snv/*.g.vcf.gz`
gvcf_commond=${gvcf_list///share/-v /share} 
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ref_genome --algo GVCFtyper $gvcf_commond 01_cohort_SNV.vcf.gz

# Left-alignment and normalized VCF
bcftools norm -m-any --check-ref w -f $ref_genome 01_cohort_SNV.vcf.gz -Oz -o 02_cohort_SNV_norm.vcf.gz
tabix 02_cohort_SNV_norm.vcf.gz



### [2] VQSR ###
mkdir -p ${outputFolderPath}/02_SNV/02_vqsr || exit 1
cd ${outputFolderPath}/02_SNV/02_vqsr

# VQSR for SNV
$GATK_INSTALL_DIR/gatk VariantRecalibrator \
	-R $ref_genome \
    -V ${outputFolderPath}/02_SNV/01_joint_calling/02_cohort_SNV_norm.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15 $snv_hapmap \
    -resource:omni,known=false,training=true,truth=true,prior=12 $snv_omni \
    -resource:1000G,known=false,training=true,truth=false,prior=10 $snv_1kg_highconf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $variant_dnsnp \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file cohort_snps.R \
	--tranches-file cohort_snps.tranches \
    -O cohort_snps.recal

$GATK_INSTALL_DIR/gatk ApplyVQSR \
	-R $ref_genome \
    -V ${outputFolderPath}/02_SNV/01_joint_calling/02_cohort_SNV_norm.vcf.gz \
    --recal-file cohort_snps.recal \
    --tranches-file cohort_snps.tranches \
    --truth-sensitivity-filter-level $vqsr_threshold \
    --create-output-variant-index true \
    -mode SNP \
    -O 01_VQSR_snp.vcf.gz

# VQSR for INDEL
$GATK_INSTALL_DIR/gatk VariantRecalibrator \
	-R $ref_genome \
    -V 01_VQSR_snp.vcf.gz \
	-resource:mills,known=false,training=true,truth=true,prior=12 $indel_mill \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $variant_dnsnp \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
    --max-gaussians 6 \
	--rscript-file cohort_indels.R \
	--tranches-file cohort_indels.tranches \
    -O cohort_indels.recal

$GATK_INSTALL_DIR/gatk ApplyVQSR \
	-R $ref_genome \
    -V 01_VQSR_snp.vcf.gz \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level $vqsr_threshold \
    --create-output-variant-index true \
    -mode INDEL \
    -O 02_VQSR_snp_indel.vcf.gz

# Filter variants labeled PASS
bcftools view -f PASS 02_VQSR_snp_indel.vcf.gz -Oz -o 03_cohort_SNV_HC.vcf.gz
tabix 03_cohort_SNV_HC.vcf.gz



### [3] Convert the big VCF to plink formats ###
mkdir -p ${outputFolderPath}/02_SNV/03_change_format || exit 1
mkdir -p ${outputFolderPath}/02_SNV/03_change_format/01_vcf_by_chr || exit 1
mkdir -p ${outputFolderPath}/02_SNV/03_change_format/02_plink || exit 1

## split vcf by chr
cd ${outputFolderPath}/02_SNV/03_change_format/01_vcf_by_chr
tabix --list-chroms ${outputFolderPath}/02_SNV/02_vqsr/03_cohort_SNV_HC.vcf.gz > chromosomes_all.txt
head -24 chromosomes_all.txt > chromosomes.txt

while IFS= read -r line; do
  tabix -h ${outputFolderPath}/02_SNV/02_vqsr/03_cohort_SNV_HC.vcf.gz $line | bgzip -c > cohort_chr$line.vcf.gz
  tabix cohort_chr$line.vcf.gz;
done < chromosomes.txt

## vcf2plink
cd ${outputFolderPath}/02_SNV/03_change_format/02_plink

$PLINK2 \
	--vcf ${outputFolderPath}/02_SNV/02_vqsr/03_cohort_SNV_HC.vcf.gz \
	--vcf-idspace-to _ \
	--double-id \
	--allow-extra-chr 0 \
	--set-missing-var-ids @:#[b37]\$r,\$a \
	--new-id-max-allele-len 1000 \
	--keep-allele-order \
	--make-bed \
	--out cohort_v1

plink --bfile cohort_v1 --allow-extra-chr 0 --split-x hg19 no-fail --keep-allele-order --make-bed --out cohort_v2
plink --bfile cohort_v2 --allow-extra-chr 0 --set-hh-missing --keep-allele-order --make-bed --out cohort_v3

# The last thing you need to change the plink.fam files manually according to the metadata file of the cohort.