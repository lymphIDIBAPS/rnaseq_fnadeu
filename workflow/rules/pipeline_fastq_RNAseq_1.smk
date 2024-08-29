#### Snakefile to do the RNAseq pipeline, provided by Ferran Nadeu ####

import argparse
import os
import subprocess
import time
import yaml



rule create_folders:
    input:
        ## Some input
        "resources/start.txt"
    output:
        ## The newly created directories
        directory(outDir)
    run:
        # Starting time
        date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")

        # Get current directory
        currentDir = os.getcwd()

        # Create output directory
        aName = "_"+options.analysisName if options.analysisName != "" else ""
        outDir = options.workDir+"/"+date_str+"_pipeline_fastq_RNAseq"+aName
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