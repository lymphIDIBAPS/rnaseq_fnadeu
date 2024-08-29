import argparse
import os
import yaml
import subprocess
import time

with open("config/config.yaml", "r") as file:
    config = yaml.safe_load(file)

date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")

# Get current directory
currentDir = os.getcwd()

# Create output directory
aName = "_" + config["arguments"]["analysisName"] if config["arguments"]["analysisName"] != "" else ""
outDir = config["arguments"]["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName
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