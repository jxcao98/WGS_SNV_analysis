#!/bin/bash

#################################################################################

## Step2: Cohort-level SNV and QC
## Usage: ./run_02_joint_calling.sh myConfigFile.config

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
mkdir -p ${outputFolderPath}/02_SNV/01_joint_calling || exit 1
runDir=${outputFolderPath}/02_SNV/01_joint_calling


# Sentieon License
export SENTIEON_LICENSE=${SENTIEON_LICENSE}

# log
logfile=${runDir}/joint_calling.log
exec >$logfile 2>&1

# =============================================================================
#                           Run Pipeline
# =============================================================================

### [1] Joint Calling ###
mkdir -p ${runDir}/01_joint_calling || exit 1
cd ${runDir}/01_joint_calling

# Joint calling
gvcf_list=`ls ${outputFolderPath}/01_sentieon/*/*/06_call_snv/*.g.vcf.gz`
gvcf_commond=${gvcf_list///share/-v /share} 
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ref_genome --algo GVCFtyper $gvcf_commond cohort_SNV.vcf.gz

# Left-alignment and normalized VCF
bcftools norm -m-any --check-ref w -f $ref_genome cohort_SNV.vcf.gz -Oz -o cohort_SNV_norm.vcf.gz
tabix cohort_SNV_norm.vcf.gz



### [2] VQSR ###
mkdir -p ${runDir}/02_vqsr || exit 1
cd ${runDir}/02_vqsr

# VQSR for SNV
$GATK_INSTALL_DIR/gatk VariantRecalibrator \
	-R $ref_genome \
    -V ${runDir}/01_joint_calling/cohort_SNV_norm.vcf.gz \
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
    -V ${runDir}/01_joint_calling/cohort_SNV_norm.vcf.gz \
    --recal-file cohort_snps.recal \
    --tranches-file cohort_snps.tranches \
    --truth-sensitivity-filter-level $vqsr_threshold \
    --create-output-variant-index true \
    -mode SNP \
    -O cohort_SNV_vqsr_snp.vcf.gz

# VQSR for INDEL
$GATK_INSTALL_DIR/gatk VariantRecalibrator \
	-R $ref_genome \
    -V cohort_SNV_vqsr_snp.vcf.gz \
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
    -V cohort_SNV_vqsr_snp.vcf.gz \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level $vqsr_threshold \
    --create-output-variant-index true \
    -mode INDEL \
    -O cohort_SNV_vqsr_all.vcf.gz

# Filter variants labeled PASS
bcftools view -f PASS cohort_SNV_vqsr_all.vcf.gz -Oz -o cohort_SNV_HC.vcf.gz
tabix cohort_SNV_HC.vcf.gz



### [3] Convert the VCF to plink formats ###
mkdir -p ${runDir}/03_vcf2plink || exit 1
cd ${runDir}/03_vcf2plink

$PLINK2 \
	--vcf ${runDir}/02_vqsr/cohort_SNV_HC.vcf.gz \
	--vcf-idspace-to _ \
	--double-id \
	--allow-extra-chr 0 \
	--set-missing-var-ids @:#[b37]\$r,\$a \
	--new-id-max-allele-len 1000 \
	--keep-allele-order \
	--make-bed \
	--out cohort_SNV_HC_v1

plink --bfile cohort_SNV_HC_v1 --allow-extra-chr 0 --split-x hg19 no-fail --keep-allele-order --make-bed --out cohort_SNV_HC_v2
plink --bfile cohort_SNV_HC_v2 --allow-extra-chr 0 --set-hh-missing --keep-allele-order --make-bed --out cohort_SNV_HC_v3

# The last thing you need to change the cohort_SNV_HC_v3.fam file manually according to the metadata file of the cohort.