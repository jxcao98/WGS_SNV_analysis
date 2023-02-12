#!/bin/bash

#################################################################################

## Step3: Quality control for cohort
## Usage: ./run_03_QC.sh myConfigFile.config

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
mkdir -p ${outputFolderPath}/02_SNV/02_QC || exit 1
runDir=${outputFolderPath}/02_SNV/02_QC

# log
logfile=${runDir}/QC.log
exec >$logfile 2>&1

# =============================================================================
#                           Run QC Pipeline
# =============================================================================
cd ${runDir}

cohort_SNV_HC=${outputFolderPath}/02_SNV/01_joint_calling/03_vcf2plink/cohort_SNV_HC_v3

### [1] Missingness ###
plink --bfile $cohort_SNV_HC --geno 0.05 --make-bed --out 01_Missingness_1
plink --bfile 01_Missingness_1 --mind 0.05 --make-bed --out 01_Missingness_2


### [2] Sex discrepancy ###
plink --bfile 01_Missingness_2 --check-sex 
Rscript --no-save $Rscript_path/gender_check.R

# Delete individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}'> sex_discrepancy.txt
# cat plink.sexcheck | awk '{if($4==1)  print}' | awk '{print$1,$2}' > sex_discrepancy.txt
plink --bfile 01_Missingness_2 --remove sex_discrepancy.txt --make-bed --out 02_Sex_coincident


### [3] Filter Variants ###
# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' 02_Sex_coincident.bim > snp_1_22.txt
plink --bfile 02_Sex_coincident --extract snp_1_22.txt --make-bed --out 03_1_autosomal_only

# Remove SNPs with a low MAF frequency.
plink --bfile 03_1_autosomal_only --maf 0.05 --make-bed --out 03_2_common_variants_only

# Remove SNPs extremely from HWE
plink --bfile 03_2_common_variants_only --hwe 1e-10 --make-bed --out 03_3_HWE_only

# LD prune
plink --bfile 03_3_HWE_only --exclude $Rscript_path/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP


### [4] Heterozygosity rate ###
plink --bfile 03_3_HWE_only --extract indepSNP.prune.in --het --out R_check
Rscript --no-save $Rscript_path/heterozygosity_outliers_list.R

# Remove heterozygosity rate outliers.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile 03_3_HWE_only --remove het_fail_ind.txt --make-bed --out 04_heterozygosity


### [5] IBD ###
# Check for relationships between individuals
$KING -b 04_heterozygosity.bed --related


