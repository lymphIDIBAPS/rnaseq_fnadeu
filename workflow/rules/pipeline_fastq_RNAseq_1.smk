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
date_str = "2024_25101412"

# Set name for all rules
aName = "_" + config["analysisName"] if config["analysisName"] != "" else ""

# Set outDir for all rules
outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName

# Set path to reference Data and pipelines
pathToReferenceDataAndPipelines = config["pathToReferenceDataAndPipelines"]
pathToScripts = "workflow/scripts/"

## Now Ferran does the scripting for the cluster, I will do it in a config file for the cluster
## Prepare and submit job script for each sample
## The command sent to each job is the following:
##
## comm = "python "+options.apathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_2.py 
## -p "+options.apathToReferenceDataAndPipelines+" -o "+outDir+" -s "+sample+" -F1 "+fq1+" -F2 "+fq2+" 
## -TS "+strand+" -tT "+options.tTrimming+" -ad "+options.adapters+" -r "+options.removeFastqs+" -rb "+options.removeBams+" 
## -sortmernaDB "+options.sortmernaDB+" -g "+options.agenes+" -qc "+options.runQC+" -@ "+options.cpus

rule create_folders:
    output:
        ## The newly created directories
        "resources/create_folders.txt"
    params:
        date_str = date_str,
        aName = aName,
    conda:
        "../envs/bash.yaml"
    shell:
        """
        bash workflow/scripts/create_folders.sh {params.date_str} {params.aName}
        """


##
## Summary QC once all previous jobs are done
## 

rule summary_qc:
    input:
        "resources/create_folders.txt",
        sortmerna_result = "{outDir}/MULTIQC_FASTQC/files/{sample}_sortmerna_1_fastqc.html",
    output:
        summary_out = "{outDir}/MULTIQC_FASTQC/files/{sample}_MULTIQC_end.txt"
    params:
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
    envmodules:
        "intel/2018.3",
        "python/3.6.5"
    shell:
        """
        touch {params.outDir}/MULTIQC_FASTQC/files/{wildcards.sample}_MULTIQC_end.txt
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
        summary_out = "{outDir}/MULTIQC_FASTQC/files/{sample}_MULTIQC_end.txt"
    output:
        plots_out = "{outDir}/MULTIQC/{sample}_pipeline_fastq_RNAseq_PCAs.pdf"
    envmodules:
        "R/3.5.1"
    params:
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
    conda:
        "../envs/plots.yaml"
    shell:
        """
        echo Lines 175 and 176 from RNAseq_1
        touch {params.outDir}/MULTIQC/{wildcards.sample}_pipeline_fastq_RNAseq_PCAs.pdf
        echo Rscript --vanilla {pathToScripts}pipeline_fastq_RNAseq_4.R {ensemblTable} {sampleTableToOpen} {kallistoPath} {outDir}/MULTIQC/{wildcards.sample}_pipeline_fastq_RNAseq_PCAs.pdf
        """

# "Rscript --vanilla "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_3.R "+ensemblTable+" "+options.infoRun+" "+outDir+"/KALLISTO "+outDir+"/MULTIQC/"+date_str+"_pipeline_fastq_RNAseq_PCAs.pdf"
