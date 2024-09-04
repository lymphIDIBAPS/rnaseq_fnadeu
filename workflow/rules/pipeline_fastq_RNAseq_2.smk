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

# sample = options.sample
# fastq1 = options.fq1
# fastq2 = options.fq2
# TranscriptionStrand = options.TranscriptionStrand
# cpus = options.cpus

# Extract sample names and paths to reads
samples = pep.sample_table

rule sortmerna:
    input:
        config_file = "config/project_config.yaml"
    output:
        "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_SORTMERNA/{sample}_aligned.log",
    params:
        forw = lambda wc: samples.loc[wc.sample]["forward"],
        reve = lambda wc: samples.loc[wc.sample]["reverse"],
        sortmerna_db = config["pathToReferenceDataAndPipelines"]+"/smr_v4.3_"+config["sortmernaDB"]+"_db.fasta",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        # outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
        cpus = config["cpus"]
    shell:
        """
        mkdir -p {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}
        sortmerna --ref {params.sortmerna_db} --reads {params.forw} --reads {params.reve} --other {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out\
        --workdir {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample} --fastx --paired_in -threads {params.cpus} -out2 -v
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned_fwd.fq.gz {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_1.fq.gz
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned_rev.fq.gz {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_2.fq.gz
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned.log {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_aligned.log
        rm -rf {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}
        """


rule multiqc_sortmerna_log:
    input:
        aligned_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_SORTMERNA/{sample}_aligned.log",
    output:
        multiqc_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/MULTIQC/files/{sample}_aligned.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        fastq1 = lambda wc: samples.loc[wc.sample]["forward"],
        fastq2 = lambda wc: samples.loc[wc.sample]["reverse"],
    run:
        # Define the paths
        logInput_path = input.aligned_log
        logTrimmomatic_path = output.multiqc_log

        # Open input and output log files
        with open(logInput_path, "r") as logInput, open(logTrimmomatic_path, "w") as logSortmerna:
            # Replace file names with sample name in the log file
            for lLine in logInput:
                logSortmerna.write(
                    lLine.replace(params.fastq1.split("/")[-1], params.sample + ".fq.gz")
                         .replace(params.fastq2.split("/")[-1], params.sample + ".fq.gz")
                )


def check_trimming(tTrimming):
    return "HEADCROP:1" if tTrimming == "yes" else ""

rule trimmomatic:
    input:
        config_file = "config/project_config.yaml"
    output:
        log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}-trimmomatic.log"
    params:
        fastq1=lambda wildcards: samples.loc[wildcards.sample]["forward"],
        fastq2=lambda wildcards: samples.loc[wildcards.sample]["reverse"],
        paired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        unpaired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_unpaired.fastq.gz",
        paired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        unpaired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_unpaired.fastq.gz",
        cpus = config["cpus"],
        adapters = "resources/TruSeq3-PE.fa",  # Default adapter file
        result = check_trimming(config["tTrimming"])   # Default trimming option
    shell:
        """
        trimmomatic PE -threads {params.cpus} -phred33 {params.fastq1} {params.fastq2} \
        {params.paired1} {params.unpaired1} {params.paired2} {params.unpaired2} \
        ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 {params.result} 2> {output.log}
        """

# java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar

rule multiqc_trimmomatic_log:
    input:
        aligned_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}-trimmomatic.log",
    output:
        multiqc_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/MULTIQC/files/{sample}-trimmomatic.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        fastq1 = lambda wc: samples.loc[wc.sample]["forward"],
        fastq2 = lambda wc: samples.loc[wc.sample]["reverse"],
    run:
        # Define the paths
        logInput_path = input.aligned_log
        logTrimmomatic_path = output.multiqc_log

        # Open input and output log files
        with open(logInput_path, "r") as logInput, open(logTrimmomatic_path, "w") as logTrimmomatic:
            # Replace file names with sample name in the log file
            for lLine in logInput:
                logTrimmomatic.write(lLine.replace("_sortmerna_1", "").replace("_sortmerna_2", ""))


rule fastqc:
    input: 
        aligned_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}-trimmomatic.log",
    output: 
        "{outDir}/MULTIQC_FASTQC/files/{sample}_sortmerna_1_fastqc.html"
    params:
        fastq1 = lambda wildcards: samples.loc[wildcards.sample]["forward"],
        fastq2 = lambda wildcards: samples.loc[wildcards.sample]["reverse"],
        fastq1_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_1.fq.gz",
        fastq2_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_2.fq.gz",
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        cpus = config["cpus"],
    shell:
        """
        fastqc -o {params.outDir}/MULTIQC_FASTQC/files/ -t {params.cpus} {params.fastq1} {params.fastq2} \
        {params.fastq1_sortmerna} {params.fastq2_sortmerna} {params.pairedFile1} {params.pairedFile2}
        """

def get_strand_flag(transcription_strand):
    if transcription_strand == "first":
        return "--fr-stranded"
    elif transcription_strand == "second":
        return "--rf-stranded"
    else:
        return ""

## kallisto_index: builds an index
## from a FASTA formatted file of target sequences. Compute intensive rule

rule kallisto_index:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        index_out_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"
    params:
        threads = int(config["cpus"]),
    shell:
        """
        kallisto index -i {output.index_out_path} --threads=24 {input.index_path}
        """


rule kallisto:
    # This rule does not work with the paired-end file samples
    input:
        "resources/start.txt",
        "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"
    output: 
        "{outDir}/KALLISTO/{sample}_abundance.h5",
        log = "{outDir}/KALLISTO/{sample}-kallisto.log"
    params:
        fastq1 = lambda wildcards: samples.loc[wildcards.sample]["forward"],
        transcription_strand = lambda wildcards: samples.loc[wildcards.sample]["TranscriptionStrand"],
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        cpus = config["cpus"],
        kallisto_ref = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx",
    run:
        # Get the strand flag based on the transcription strand
        # From the kallisto quant I deleted the {strand_flag} because it did not find any pseudoalignment with it
        strand_flag = get_strand_flag(params.transcription_strand)
        
        shell ("""
        mkdir -p {params.outDir}/KALLISTO/{wildcards.sample}
        touch {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5
        kallisto quant -t {params.cpus} -i {params.kallisto_ref} -o {params.outDir}/KALLISTO/{wildcards.sample} {params.pairedFile1} {params.pairedFile2} 2> {output.log}
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5 {params.outDir}/KALLISTO/{wildcards.sample}_abundance.h5
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.tsv {params.outDir}/KALLISTO/{wildcards.sample}_abundance.tsv
        mv {params.outDir}/KALLISTO/{wildcards.sample}/run_info.json {params.outDir}/KALLISTO/{wildcards.sample}_run_info.json
        """)

rule multiqc_kallisto_log:
    input:
        kallisto_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/KALLISTO/{sample}-kallisto.log"
    output:
        multiqc_log = "{outDir}/MULTIQC/files/{sample}-kallisto.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        fastq1 = lambda wc: samples.loc[wc.sample]["forward"],
        fastq2 = lambda wc: samples.loc[wc.sample]["reverse"],
    run:
        # Define the paths
        logInput_path = input.kallisto_log
        logKallisto_path = output.multiqc_log

        # Open input and output log files
        with open(logInput_path, "r") as logInput, open(logKallisto_path, "w") as logKallisto:
            # Replace file names with sample name in the log file
            for lLine in logInput:
                logKallisto.write(lLine.replace("_sortmerna_1", "").replace("_sortmerna_2", ""))