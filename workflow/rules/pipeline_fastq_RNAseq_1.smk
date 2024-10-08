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
pathToReferenceDataAndPipelines = config["pathToReferenceDataAndPipelines"]
pathToScripts = "workflow/scripts/"

rule create_folders:
    input:
        ## Some input
        "resources/start.txt"
    output:
        ## The newly created directories
        "resources/create_folders.txt"
    params:
        date_str = date_str,
        aName = aName,
    shell:
        """
        bash workflow/scripts/create_folders.sh {params.date_str} {params.aName}
        """


##
## Summary QC once all previous jobs are done
## 

rule summary_qc:
    input:
        "resources/create_folders.txt"
    output:
        f"resources/{date_str}_MULTIQC_end.txt"
    envmodules:
        "python/3.6.5"
        "gcc/9.2.0"
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

ensemblTable = pathToReferenceDataAndPipelines + "/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38_genes_transcripts_release-105.txt"

sampleTableToOpen = "config/project_config.yaml"

kallistoPath = outDir+"/KALLISTO"

rule plots:
    input:
        f"{outDir}/MULTIQC/{date_str}_pipeline_fastq_RNAseq_BAM{aName}_data/kallisto_alignment.txt"
    output:
        f"{outDir}'/MULTIQC/'{date_str}'_pipeline_fastq_RNAseq_PCAs.pdf"
        # f"/resources/{date_str}_plots_made.txt"
    envmodules:
        "python/3.6.5"
        "zlib/1.2.11" 
        "szip/2.1.1"
        "R/3.5.0"
    shell:
        """
        echo Lines 175 and 176 from RNAseq_1
        touch {outDir}/MULTIQC/{date_str}_pipeline_fastq_RNAseq_PCAs.pdf
        echo Rscript --vanilla {pathToScripts}pipeline_fastq_RNAseq_4.R {ensemblTable} {sampleTableToOpen} {kallistoPath} {outDir}/MULTIQC/{date_str}_pipeline_fastq_RNAseq_PCAs.pdf
        """

# "Rscript --vanilla "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_3.R "+ensemblTable+" "+options.infoRun+" "+outDir+"/KALLISTO "+outDir+"/MULTIQC/"+date_str+"_pipeline_fastq_RNAseq_PCAs.pdf"
