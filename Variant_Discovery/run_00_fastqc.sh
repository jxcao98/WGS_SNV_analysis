#!/bin/bash

#################################################################################

## Step0: QC for original fastq files
## Usage: ./run_00_fastqc.sh myConfigFile.config

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
mkdir -p ${outputFolderPath}/00_fastqc || exit 1
mkdir -p ${outputFolderPath}/00_fastqc/${batch} || exit 1
mkdir -p ${outputFolderPath}/00_fastqc/${batch}/log || exit 1
mkdir -p ${outputFolderPath}/00_fastqc/${batch}/tmp || exit 1
# mkdir -p ${outputFolderPath}/${batch}/log || exit 1


# =============================================================================
#                           Run QC Pipeline
# =============================================================================
## Loop over each sample
while [ 1 ]
do

	# Sample name and Fastq path
	read baminfo || break
	array=(${baminfo//,/ })
	sample_name=${array[0]}
	fastq_1=${array[1]}
	fastq_2=${array[2]}
		
	# log
	logfile=${outputFolderPath}/00_fastqc/${batch}/log/${sample_name}.log
	exec >$logfile 2>&1
	
	# run fastqc
	$fastqc -t $nt $fastq_1 -o ${outputFolderPath}/00_fastqc/${batch}/ -d ${outputFolderPath}/00_fastqc/${batch}/tmp
	$fastqc -t $nt $fastq_2 -o ${outputFolderPath}/00_fastqc/${batch}/ -d ${outputFolderPath}/00_fastqc/${batch}/tmp
	
done < $id_mapping

# run multiqc
cd ${outputFolderPath}/00_fastqc/${batch}/
$multiqc ${outputFolderPath}/00_fastqc/${batch}/
