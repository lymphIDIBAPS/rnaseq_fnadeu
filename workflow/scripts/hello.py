import os
import yaml
# import subprocess
import time

with open("config/config.yaml", "r") as file:
    config = yaml.safe_load(file)

date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")

# Get current directory
currentDir = os.getcwd()

# Create output directory
aName = "_" + config["analysisName"] if config["analysisName"] != "" else ""
outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName
os.makedirs(outDir)

# Create subfolders
os.makedirs(outDir+"/log")
os.makedirs(outDir+"/FASTQ")
os.makedirs(outDir+"/FASTQ_SORTMERNA")
os.makedirs(outDir+"/FASTQ_TRIMMED")
os.makedirs(outDir+"/BAM")
os.makedirs(outDir+"/KALLISTO")
os.makedirs(outDir+"/MULTIQC_FASTQC")
os.makedirs(outDir+"/MULTIQC_FASTQC/files")
os.makedirs(outDir+"/MULTIQC")
os.makedirs(outDir+"/MULTIQC/files")


import peppy



import peppy

# Load the PEP file
pep = peppy.Project("config/project_config.yaml")

# Extract sample names and paths to reads
samples = pep.sample_table

# Extract samples paths to reads
samples2 = pep.sample_table["sample_name"]

# Function to look up forward and reverse reads
def get_fastq_paths(wildcards, strand):
    sample_info = samples.loc[wildcards.sample]
    return sample_info[strand]

