# WGS_SNV_analysis
## WGS variant discovery and quality control
Before executing this pipeline, you need to create a comma-separated file to specify the sample name and fastq file path (Sample ID, Fastq1, Fastq2). Headers is not needed.
### Step0 Check the quality of the sequence data
`sh run_00_fastqc.sh myConfigFile.config`
### Step1 Perform sequence alignment and SNV detection on individual samples
In this step, sentieon is used to perform sequence alignment and variant detection, which is a commercial software used as an alternative to the bwa+gatk process.  
  
  
`sh run_01_sentieon.sh myConfigFile.config`
### Step2 Perform cohort-level variant detection
`sh run_02_joint_calling.sh myConfigFile.config`
### Step3 Quality control to exclude low quality variants and samples
`sh run_03_QC.sh myConfigFile.config`
### Step4 Variant annotations for prioritization
`sh run_04_Annotation.sh myConfigFile.config`
