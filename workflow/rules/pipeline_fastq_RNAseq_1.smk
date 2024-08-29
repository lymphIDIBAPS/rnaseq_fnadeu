#### Snakefile to do the RNAseq pipeline, provided by Ferran Nadeu ####

import argparse
import os
import subprocess
import time
import yaml

# Open configfile to access arguments:
configfile: "config/config.yaml"

rule create_folders:
    input:
        ## Some input
        "resources/start.txt"
    output:
        ## The newly created directories
        # directory(outDir)
        "resources/end.txt"
    run:
        # Open the config file to load options
        with open("config/config.yaml", "r") as file:
            config = yaml.safe_load(file)
        
        # Starting time
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

        # Create end.txt
        fp = open("resources/end.txt", 'w')
        fp.close()


## Now Ferran does the scripting for the cluster, I will do it in a config file for the cluster
## Prepare and submit job script for each sample
## The command sent to each job is the following:
##
## comm = "python "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_2.py 
## -p "+options.pathToReferenceDataAndPipelines+" -o "+outDir+" -s "+sample+" -F1 "+fq1+" -F2 "+fq2+" 
## -TS "+strand+" -tT "+options.tTrimming+" -ad "+options.adapters+" -r "+options.removeFastqs+" -rb "+options.removeBams+" 
## -sortmernaDB "+options.sortmernaDB+" -g "+options.genes+" -qc "+options.runQC+" -@ "+options.cpus

##
## Summary QC and plots once all previous jobs are done
## git commit -m "Add rule summary_qc"

rule summary_qc:
    input:
        "resources/end.txt"
    output:
        "resources/MULTIQC_end.txt"
    run:
        date_str

        # Create end.txt
        fp = open("resources/MULTIQC_end.txt", 'w')
        fp.close()