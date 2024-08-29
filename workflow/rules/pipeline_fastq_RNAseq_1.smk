#### Snakefile to do the RNAseq pipeline, provided by Ferran Nadeu ####

import argparse
import os
import subprocess
import time
import yaml

# Open configfile to access arguments:
configfile: "config/config.yaml"

# Set date for all rules
date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")

# Set name for all rules
aName = "_" + config["analysisName"] if config["analysisName"] != "" else ""

# Set outDir for all rules
outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName

rule create_folders:
    input:
        ## Some input
        "resources/start.txt"
    output:
        ## The newly created directories
        "resources/create_folders.txt"
    run:
        # Open the config file to load options
        with open("config/config.yaml", "r") as file:
            config = yaml.safe_load(file)

        # Get current directory
        currentDir = os.getcwd()

        # Create output directory
        # outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName
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

        # Make create_folders.txt
        fp = open("resources/create_folders.txt", 'w')
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
## Summary QC once all previous jobs are done
## 

rule summary_qc:
    input:
        "resources/create_folders.txt"
    output:
        f"resources/{date_str}_MULTIQC_end.txt"
    shell:
        """
        echo {date_str}
        touch "resources/{date_str}_MULTIQC_end.txt"
        echo {config[workDir]}
        multiqc -f -i "${date_str}_pipeline_fastq_RNAseq_FASTQC${aName}" -b 'Multiqc report for RNAseq pipeline (FASTQC)' -n "${date_str}_pipeline_fastq_RNAseq_FASTQC${aName}" -o "${outDir}/MULTIQC_FASTQC" "${outDir}/MULTIQC_FASTQC/files"
        multiqc -f -i "${date_str}_pipeline_fastq_RNAseq_BAM${aName}" -b 'Multiqc report for RNAseq pipeline (BAM)' -n "${date_str}_pipeline_fastq_RNAseq_BAM${aName}" -o "${outDir}/MULTIQC" "${outDir}/MULTIQC/files"
        """

##
## Summary plots once all previous jobs are done
## 

rule plots:
    input:
        "/resources/create_folders.txt"
    output:
        f"{outDir}'/MULTIQC/'{date_str}'_pipeline_fastq_RNAseq_PCAs.pdf"
    shell:
        """
        # Lines 175 and 176 from RNAseq_2
        """