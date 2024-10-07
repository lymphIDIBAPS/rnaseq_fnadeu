#!/bin/bash

# Load config.yaml
config_file="config/config.yaml"

# Extract workDir, remove quotes
workDir=$(grep "workDir" $config_file | awk '{print $2}' | sed 's/"//g')

# Get current directory
currentDir=$(pwd)

# Input variables (date_str and aName passed as script arguments)
date_str=$1  # First argument for date string
aName=$2     # Second argument for aName

# Create output directory
outDir="${workDir}/${date_str}_pipeline_fastq_RNAseq${aName}"
mkdir -p "$outDir"

# Create subfolders
mkdir -p "$outDir/log"
mkdir -p "$outDir/FASTQ"
mkdir -p "$outDir/FASTQ_SORTMERNA"
mkdir -p "$outDir/FASTQ_TRIMMED"
mkdir -p "$outDir/BAM"
mkdir -p "$outDir/KALLISTO"
mkdir -p "$outDir/MULTIQC_FASTQC/files"
mkdir -p "$outDir/MULTIQC/files"

# Create an empty create_folders.txt file
touch resources/create_folders.txt
