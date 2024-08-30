#### Snakefile to do the RNAseq pipeline, provided by Ferran Nadeu ####

import argparse
import os
import subprocess
import time
import yaml
import peppy
import pandas as pd

# PEP file with info on our project
pepfile: "config/project_config.yaml"

# Open configfile to access arguments
configfile: "config/config.yaml"


## Now Ferran does the scripting for the cluster, I will do it in a config file for the cluster
## Prepare and submit job script for each sample
## The command sent to each job is the following:
##
## comm = "python "+options.apathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_2.py 
## -p "+options.apathToReferenceDataAndPipelines+" -o "+outDir+" -s "+sample+" -F1 "+fq1+" -F2 "+fq2+" 
## -TS "+strand+" -tT "+options.tTrimming+" -ad "+options.adapters+" -r "+options.removeFastqs+" -rb "+options.removeBams+" 
## -sortmernaDB "+options.sortmernaDB+" -g "+options.agenes+" -qc "+options.runQC+" -@ "+options.cpus


# 0) define variables and files

sortmerna = config["pathToReferenceDataAndPipelines"]+"/programs/sortmerna-4.3.4/bin/sortmerna"
sortmerna_db = config["pathToReferenceDataAndPipelines"]+"/programs/sortmerna-4.3.4/db/smr_v4.3_"+config["sortmernaDB"]+"_db.fasta"
star_genome = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/STAR_100"
refFlat = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GFF3/Homo_sapiens.GRCh38.105.refFlat"
ribosomal_intervals = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GTF/Homo_sapiens.GRCh38.105.ribosomalIntervals"
rseqc_bed = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GFF3/Homo_sapiens.GRCh38.105.bed"

if config["genes"] == "cDNA": 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.release-105.idx"
elif config["genes"] == "ncRNA": 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.ncrna.release-105.idx"
else: 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.ncrna.release-105.idx"

outDir = config["workDir"]
# sample = options.sample
# fastq1 = options.fq1
# fastq2 = options.fq2
# TranscriptionStrand = options.TranscriptionStrand
# cpus = options.cpus

rule sortmerna:
    input:  "config/project_config.yaml"
    output: "results/{sample}.txt"
    shell:
        """
        touch results/{wildcards.sample}.txt
        """