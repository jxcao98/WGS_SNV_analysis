# WGS_SNV_analysis
## WGS variant discovery and quality control
Before executing this pipeline, you need to create a comma-separated file to specify the sample name and fastq file path (Sample ID, Fastq1, Fastq2). Headers is not needed. An example is `/Variant_Discovery/cohort_info.txt`. You also need a config files to specify necessary path of tools and data, like `/Variant_Discovery/myConfigFile.config`
### Step0 Check the quality of the sequence data
`sh /Variant_Discovery/run_00_fastqc.sh myConfigFile.config`
### Step1 Perform sequence alignment and SNV detection on individual samples
In this step, sentieon is used to perform sequence alignment and variant detection, which is a commercial software used as an alternative to the bwa+gatk process.  
  
  
`sh /Variant_Discovery/run_01_sentieon.sh myConfigFile.config`
### Step2 Perform cohort-level variant detection
`sh /Variant_Discovery/run_02_joint_calling.sh myConfigFile.config`
### Step3 Quality control to exclude low quality variants and samples
`sh /Variant_Discovery/run_03_QC.sh myConfigFile.config`
### Step4 Variant annotations for prioritization
`sh /Variant_Discovery/run_04_Annotation.sh myConfigFile.config`
