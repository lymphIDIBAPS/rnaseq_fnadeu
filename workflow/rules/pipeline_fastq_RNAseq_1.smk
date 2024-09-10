#### Snakefile to do the RNAseq pipeline, provided by Ferran Nadeu ####

import argparse
import os
import subprocess
import time
import yaml

# Open configfile to access arguments:
configfile: "config/config.yaml"

# Set date for all rules
# date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")
date_str = "20240902_162657"

# Set name for all rules
aName = "_" + config["analysisName"] if config["analysisName"] != "" else ""

# Set outDir for all rules
outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName

# Set path to reference Data and pipelines
# pathToReferenceDataAndPipelines = config["pathToReferenceDataAndPipelines"]
pathToReferenceDataAndPipelines = "workflow/scripts/"

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
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName
        os.makedirs(outDir, exist_ok=True)

        # Create subfolders
        os.makedirs(outDir+"/log", exist_ok=True)
        os.makedirs(outDir+"/FASTQ", exist_ok=True)
        os.makedirs(outDir+"/FASTQ_SORTMERNA", exist_ok=True)
        os.makedirs(outDir+"/FASTQ_TRIMMED", exist_ok=True)
        os.makedirs(outDir+"/BAM", exist_ok=True)
        os.makedirs(outDir+"/KALLISTO", exist_ok=True)
        os.makedirs(outDir+"/MULTIQC_FASTQC", exist_ok=True)
        os.makedirs(outDir+"/MULTIQC_FASTQC/files", exist_ok=True)
        os.makedirs(outDir+"/MULTIQC", exist_ok=True)
        os.makedirs(outDir+"/MULTIQC/files", exist_ok=True)

        # Make create_folders.txt
        fp = open("resources/create_folders.txt", 'w')
        fp.close()


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
        multiqc -f -i "{date_str}_pipeline_fastq_RNAseq_FASTQC{aName}" -b 'Multiqc report for RNAseq pipeline (FASTQC)' -n "{date_str}_pipeline_fastq_RNAseq_FASTQC{aName}" -o "{outDir}/MULTIQC_FASTQC" "{outDir}/MULTIQC_FASTQC/files"
        multiqc -f -i "{date_str}_pipeline_fastq_RNAseq_BAM{aName}" -b 'Multiqc report for RNAseq pipeline (BAM)' -n "{date_str}_pipeline_fastq_RNAseq_BAM{aName}" -o "{outDir}/MULTIQC" "{outDir}/MULTIQC/files"
        """

##
## Summary plots once all previous jobs are done
##

ensemblTable = pathToReferenceDataAndPipelines + "genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38_genes_transcripts_release-105.txt"

rule plots:
    input:
        "resources/create_folders.txt"
    output:
        # f"{outDir}'/MULTIQC/'{date_str}'_pipeline_fastq_RNAseq_PCAs.pdf"
        f"/resources/{date_str}_plots_made.txt"
    shell:
        """
        echo Lines 175 and 176 from RNAseq_1
        echo {pathToReferenceDataAndPipelines} {ensemblTable}
        touch "resources/{date_str}_plots_made.txt"
        Rscript --vanilla {pathToReferenceDataAndPipelines}pipeline_fastq_RNAseq_3.R {ensemblTable} 
        """

# "Rscript --vanilla "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_3.R "+ensemblTable+" "+options.infoRun+" "+outDir+"/KALLISTO "+outDir+"/MULTIQC/"+date_str+"_pipeline_fastq_RNAseq_PCAs.pdf"
