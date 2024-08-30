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

# Load the PEP file
pep = peppy.Project("config/project_config.yaml")

# Extract sample names and paths to reads
samples = pep.sample_table

# Function to look up forward and reverse reads
def get_fastq_paths(wildcards, strand):
    sample_info = samples.loc[wildcards.sample]
    return sample_info[strand]

# rule sortmerna:
#     input:  "config/project_config.yaml"
#     output: "results/{sample}.txt"
#     params:
#         # fastq = samplesheet("TranscriptionStrand")
#         fastq = lambda wc: lookup_sample_table(sample = wc.sample, target = "TranscriptionStrand")
#     shell:
#         """
#         touch results/{wildcards.sample}.txt
#         echo {params.fastq}
#         """

rule sortmerna_GPT:
    input:
        config_file="config/project_config.yaml",
        forw=lambda wc: get_fastq_paths(wc, "forward"),
        reve=lambda wc: get_fastq_paths(wc, "reverse"),
    output:
        "results/{sample}.txt"
    shell:
        """
        echo {input.forw} {input.reve}
        """
    
# /slgpfs/home/cli84/cli84075/bin/sortmerna \
# --ref ./sortmernaDB/silva-bac-16s-id90.fasta \
# --ref ./sortmernaDB/silva-bac-23s-id98.fasta \
# --ref ./sortmernaDB/silva-arc-16s-id95.fasta \
# --ref ./sortmernaDB/silva-arc-23s-id98.fasta \
# --ref ./sortmernaDB/silva-euk-18s-id95.fasta \
# --ref ./sortmernaDB/silva-euk-28s-id98.fasta \
# --ref ./sortmernaDB/rfam-5s-database-id98.fasta \
# --ref ./sortmernaDB/rfam-5.8s-database-id98.fasta \
# --reads {input.forward} \
# --reads {input.reverse} \
# --aligned results/{wildcards.sample}_aligned \
# --other results/{wildcards.sample}_other \
# --workdir results/{wildcards.sample}_workdir \
# --fastx \
# --paired_in \
# --threads 16 \
# --out2 \
# -v \
# --idx-dir ./sortmernaDB/index2 > {output}