#!/bin/bash

#################################################################################

## Step1: Sequence alignment and SNV detection using sentieon
## Usage: ./run_01_sentieon.sh myConfigFile.config

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
mkdir -p ${outputFolderPath}/01_sentieon || exit 1
mkdir -p ${outputFolderPath}/01_sentieon/${batch} || exit 1
mkdir -p ${outputFolderPath}/01_sentieon/${batch}/tmp || exit 1
# mkdir -p ${outputFolderPath}/${cohort}/log || exit 1

export SENTIEON_LICENSE=${SENTIEON_LICENSE}

# =============================================================================
#                           Run Sentieon Pipeline
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
	
	# outputFolderPath fold
	sampleFolderPath=${outputFolderPath}/01_sentieon/${batch}/${sample_name}
	mkdir -p ${sampleFolderPath}
	
	# log
	logfile=${sampleFolderPath}/${sample_name}.log
	exec >$logfile 2>&1
	
	
	### [1] Map reads to reference
	mkdir $sampleFolderPath/01_alignment
	cd $sampleFolderPath/01_alignment
	
	echo $sample_name
	( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\tID:$sample_name\tSM:$sample_name\tPL:ILLUMINA" \
		-t $nt -K 10000000 $ref_genome $fastq_1 $fastq_2 || echo -n 'error' ) \
		| $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $ref_genome -o $sample_name.sorted.bam -t $nt --sam2bam -i -


	### [2] Calculate data metrics (Optional)


	### [3] Remove or mark duplicates
	mkdir $sampleFolderPath/03_mark_dup
	cd $sampleFolderPath/03_mark_dup
	
	$SENTIEON_INSTALL_DIR/bin/sentieon driver \
		-t $nt \
		-i $sampleFolderPath/01_alignment/$sample_name.sorted.bam \
		--algo LocusCollector \
		--fun score_info $sample_name.score.txt

	$SENTIEON_INSTALL_DIR/bin/sentieon driver \
		-t $nt \
		-i $sampleFolderPath/01_alignment/$sample_name.sorted.bam \
		--algo Dedup \
		--rmdup \
		--score_info $sample_name.score.txt \
		--metrics $sample_name.dedup_metrics.txt $sample_name.deduped.bam 
	
	
	### [4] Indel realignment (optional)
	
	
	### [5] Base quality score recalibration
	mkdir $sampleFolderPath/05_BQSR
	cd $sampleFolderPath/05_BQSR
	
	$SENTIEON_INSTALL_DIR/bin/sentieon driver \
		-r $ref_genome \
		-t $nt \
		-i $sampleFolderPath/03_mark_dup/$sample_name.deduped.bam \
		--algo QualCal \
		-k $dbsnp \
		-k $known_Mills_indels \
		-k $known_1000G_indels \
		$sample_name.recal_data.table

	
	### [6] Variant calling
	mkdir $sampleFolderPath/06_call_snv
	cd $sampleFolderPath/06_call_snv

	$SENTIEON_INSTALL_DIR/bin/sentieon driver \
		-r $ref_genome \
		-t $nt \
		-i $sampleFolderPath/03_mark_dup/$sample_name.deduped.bam \
		-q $sampleFolderPath/05_BQSR/$sample_name.recal_data.table \
		--algo Haplotyper \
		-d $dbsnp \
		--emit_mode gvcf \
		$sample_name.g.vcf.gz
		
done < $id_mapping
