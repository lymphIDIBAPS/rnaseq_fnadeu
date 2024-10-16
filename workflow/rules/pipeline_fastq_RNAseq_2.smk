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
fastaFiles = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_DNA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtfFile= config["pathToReferenceDataAndPipelines"] +"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GTF/Homo_sapiens.GRCh38.105.gtf"
genomeDir = config["pathToReferenceDataAndPipelines"] +"/data/genome_GRCh38.p13_GCA_000001405.28/STAR_100/"
refFlat = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GFF3/Homo_sapiens.GRCh38.105.refFlat"
ribosomal_intervals = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GTF/Homo_sapiens.GRCh38.105.ribosomalIntervals"
THREADS = int(config["cpus"])

# Adapters for Trimmomatic
if config["adapters"] == "illumina":
	adaptersSeq = config["pathToReferenceDataAndPipelines"]+"/TruSeq3-PE-2.fa"
elif config["adapters"] == "bioskryb":
	adaptersSeq = config["pathToReferenceDataAndPipelines"]+"/Bioskryb_ResolveOME.fa"

# Index file for kallisto including only cDNA, ncRNA, or both
if config["genes"] == "cDNA":
    kallisto_path = config["pathToReferenceDataAndPipelines"] +"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_cDNA/Homo_sapiens.GRCh38.cdna.all.fa"
    kallisto_index = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.release-105.idx"
elif config["genes"] == "ncRNA":
    kallisto_path = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_ncRNA/Homo_sapiens.GRCh38.ncrna.fa"
    kallisto_index = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.ncrna.release-105.idx"
else:
    kallisto_index = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.ncrna.release-105.idx"

# Extract sample names and paths to reads
samples = pep.sample_table

### Functions ###

def check_trimming(tTrimming):
    return "HEADCROP:1" if tTrimming == "yes" else ""

# Get the strand flag for rule kallisto based on the transcription strand
def get_strand_flag(transcription_strand):
    if transcription_strand == "first":
        return "--fr-stranded"
    elif transcription_strand == "second":
        return "--rf-stranded"
    else:
        return ""

def get_transcription_strand(transcription_strand):
    if transcription_strand == "first":
        return "FIRST_READ_TRANSCRIPTION_STRAND"
    elif transcription_strand == "second":
        return "SECOND_READ_TRANSCRIPTION_STRAND"
    else:
        return "NONE"


#### PIPELINE RULES ####

rule sortmerna:
    input:
        config_file = "config/project_config.yaml"
    output:
        "{outDir}/FASTQ_SORTMERNA/{sample}_aligned.log",
    params:
        forw = lambda wc: samples.loc[wc.sample]["forward"],
        reve = lambda wc: samples.loc[wc.sample]["reverse"],
        sortmerna_db = config["pathToReferenceDataAndPipelines"]+"/programs/sortmerna-4.3.4/db/smr_v4.3_"+config["sortmernaDB"]+"_db.fasta",
        # outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
    threads:
        THREADS
    envmodules:
        "sortmerna/4.3.4"
    conda:
        "../envs/sortmerna.yaml"
    shell:
        """
        mkdir -p {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}
        sortmerna --ref {params.sortmerna_db} --reads {params.forw} --reads {params.reve} --other {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out\
        --workdir {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample} --fastx --paired_in -threads {threads} -out2 -v 
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned_fwd.fq.gz {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_1.fq.gz
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned_rev.fq.gz {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_2.fq.gz
        mv {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}/out/aligned.log {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}_aligned.log
        rm -rf {params.outDir}/FASTQ_SORTMERNA/{wildcards.sample}
        """


rule multiqc_sortmerna_log:
    input:
        aligned_log = "{outDir}/FASTQ_SORTMERNA/{sample}_aligned.log",
    output:
        multiqc_log = "{outDir}/MULTIQC/files/{sample}_aligned.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
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


rule trimmomatic:
    input:
        config_file = "config/project_config.yaml"
    output:
        log = "{outDir}/FASTQ_TRIMMED/{sample}_trimmomatic.log"
    params:
        fastq1=lambda wildcards: samples.loc[wildcards.sample]["forward"],
        fastq2=lambda wildcards: samples.loc[wildcards.sample]["reverse"],
        paired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        unpaired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_unpaired.fastq.gz",
        paired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        unpaired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_unpaired.fastq.gz",
        adapters = adaptersSeq,
        result = check_trimming(config["tTrimming"])
    threads:
        THREADS
    envmodules:
        "java/12.0.2"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar PE -threads {threads} -phred33 {params.fastq1} {params.fastq2} \
        {params.paired1} {params.unpaired1} {params.paired2} {params.unpaired2} \
        ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 {params.result} 2> {output.log}
        """
# If we want to run it in a local machine, we need to substitute "java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar"
# by only "trimmomatic"


rule multiqc_trimmomatic_log:
    input:
        aligned_log = "{outDir}/FASTQ_TRIMMED/{sample}_trimmomatic.log",
    output:
        multiqc_log = "{outDir}/MULTIQC/files/{sample}_trimmomatic.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
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
        aligned_log = "{outDir}/FASTQ_TRIMMED/{sample}_trimmomatic.log",
    output: 
        "{outDir}/MULTIQC_FASTQC/files/{sample}_sortmerna_1_fastqc.html"
    params:
        fastq1 = lambda wildcards: samples.loc[wildcards.sample]["forward"],
        fastq2 = lambda wildcards: samples.loc[wildcards.sample]["reverse"],
        fastq1_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_1.fq.gz",
        fastq2_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_2.fq.gz",
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
    threads:
        THREADS
    envmodules:
        "fastqc/0.11.9"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc -o {params.outDir}/MULTIQC_FASTQC/files/ -t {threads} {params.fastq1} {params.fastq2} \
        {params.fastq1_sortmerna} {params.fastq2_sortmerna} {params.pairedFile1} {params.pairedFile2}
        """


## kallisto_index: builds an index
## from a FASTA formatted file of target sequences. Compute intensive rule

rule kallisto_index:
    input:
        kallisto_gen = kallisto_path
    output:
        index_out_path = kallisto_index
    threads:
        THREADS
    envmodules:
        "intel/2018.3",
        "impi/2018.3",
        "zlib/1.2.11",
        "gcc/12.2.0",
        "hdf5/1.10.2",
        "szip/2.1.1",
        "kallisto/0.46.1"
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index -i {output.index_out_path} --threads={threads} {input.kallisto_gen}
        """


rule kallisto:
    # This rule does not work with the paired-end file samples
    input:
        "resources/start.txt",
        kallisto_index
    output: 
        "{outDir}/KALLISTO/{sample}_abundance.h5",
        log = "{outDir}/KALLISTO/{sample}_kallisto.log"
    params:
        fastq1 = lambda wildcards: samples.loc[wildcards.sample]["forward"],
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
        transcription_strand = get_strand_flag(config["TranscriptionStrand"]),
        kallisto_index = kallisto_index,
    threads:
        THREADS
    envmodules:
        "intel/2018.3",
        "impi/2018.3",
        "zlib/1.2.11",
        "gcc/12.2.0",
        "hdf5/1.10.2",
        "szip/2.1.1",
        "kallisto/0.46.1"
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        mkdir -p {params.outDir}/KALLISTO/{wildcards.sample}
        touch {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5
        kallisto quant -t {threads} -i {params.kallisto_index} -o {params.outDir}/KALLISTO/{wildcards.sample} {params.pairedFile1} {params.pairedFile2} {params.transcription_strand}\
        2> {output.log}
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5 {params.outDir}/KALLISTO/{wildcards.sample}_abundance.h5
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.tsv {params.outDir}/KALLISTO/{wildcards.sample}_abundance.tsv
        mv {params.outDir}/KALLISTO/{wildcards.sample}/run_info.json {params.outDir}/KALLISTO/{wildcards.sample}_run_info.json
        rm -rf {params.outDir}/KALLISTO/{wildcards.sample}
        """        


rule multiqc_kallisto_log:
    input:
        kallisto_log = "{outDir}/KALLISTO/{sample}_kallisto.log"
    output:
        multiqc_log = "{outDir}/MULTIQC/files/{sample}_kallisto.log"
    params:
        sample = "{sample}",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
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


### From here, rules are conditional on runQC: yes

rule star_genome:
    output: 
        "{genomeDir}/transcriptInfo.tab"
    params:
        fastaFiles = fastaFiles,
        gtfFile = gtfFile,
        genomeDir = genomeDir,
    threads:
        THREADS
    envmodules:
        "intel/2018.3",
        "star/2.7.8a"
    conda:
        "../envs/star.yaml"
    shell:
        """
        mkdir -p {params.genomeDir}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genomeDir} --genomeFastaFiles {params.fastaFiles} --sjdbGTFfile {params.gtfFile} --sjdbOverhang 100
        """


rule star_map:
    input:
        f"{genomeDir}transcriptInfo.tab"
    output: 
        "{outDir}/BAM/{sample}_Aligned.out.bam"
    params:
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
        genomeDir = genomeDir,
    threads:
        THREADS
    envmodules:
        "intel/2018.3",
        "star/2.7.8a"
    conda:
        "../envs/star.yaml"
    # log:
    #     "{outDir}/BAM/{sample}_Log.out"
    shell:
        """
        mkdir -p {params.genomeDir}
        STAR --runThreadN {threads} --genomeDir {params.genomeDir} --readFilesIn {params.pairedFile1} {params.pairedFile2} \
        --readFilesCommand zcat --outFileNamePrefix {params.outDir}/BAM/{wildcards.sample}_ --outSAMtype BAM Unsorted --twopassMode Basic --outBAMcompression 10
        """
# not used now compared to Romina's code: --outSAMstrandField intronMotif


rule samtools:
    input:
        unsorted_bam = "{outDir}/BAM/{sample}_Aligned.out.bam",
    output:
        sorted_bam = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam",
        sorted_bam_bai = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam.bai",
    params:
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
        star_genome = "{outDir}/BAM/{sample}__STARgenome/",
        star_pass1 = "{outDir}/BAM/{sample}__STARpass1/",
        log_final_out = "{outDir}/BAM/{sample}_Log.final.out",
        log_final_multiqc = "{outDir}/MULTIQC/files/{sample}_Log.final.out",
    threads:
        THREADS
    envmodules:
        "samtools/1.9"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -m 3G -o {output.sorted_bam} {input.unsorted_bam}
        samtools index -@ {threads} {output.sorted_bam} -o {output.sorted_bam_bai}
        rm -rf {params.star_genome} {params.star_pass1}
        cp {params.log_final_out} {params.log_final_multiqc}
        """


rule generate_md5sum:
    input:
        bam = "{outDir}/BAM/{sample}_Aligned.out.bam",
        sorted_bam = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam",
        sorted_bam_bai = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam.bai",
    output:
        bam_md5 = "{outDir}/BAM/{sample}_Aligned.out.bam.md5",
        sorted_bam_md5 = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam.md5",
        sorted_bam_bai_md5 = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam.bai.md5",
    shell:
        """
        md5sum {input.bam} > {output.bam_md5}
        md5sum {input.sorted_bam} > {output.sorted_bam_md5}
        md5sum {input.sorted_bam_bai} > {output.sorted_bam_bai_md5}
        """


rule collectRNASeqMetrics:
    # This rule fails due to not having the corresponding files
    output:
        "{outDir}/MULTIQC/files/{sample}.CollectRnaSeqMetrics"
    params:
        outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName,
        transcription_strand = get_transcription_strand(config["TranscriptionStrand"]),
        refFlat = refFlat,
        ribosomal_intervals = ribosomal_intervals
    envmodules:
        "java/12.0.2",
        "picard/2.24.0"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectRnaSeqMetrics -I {params.outDir}/BAM/{wildcards.sample}_Aligned.out.sorted.bam \
        -O {params.outDir}/MULTIQC/files/{wildcards.sample}.CollectRnaSeqMetrics -REF_FLAT {params.refFlat} -STRAND {params.transcription_strand}\
         -RIBOSOMAL_INTERVALS {params.ribosomal_intervals}
        """


rule samtools_multiqc:
    input:
        sorted_bam = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam",
    output:
        idxstats = "{outDir}/MULTIQC/files/{sample}.idxstats",
        flagstat = "{outDir}/MULTIQC/files/{sample}.flagstat",
        stats = "{outDir}/MULTIQC/files/{sample}.stats",
    threads:
        THREADS
    envmodules:
        "samtools/1.9"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats -@ {threads} {input.sorted_bam} > {output.idxstats}
        samtools flagstat -@ {threads} {input.sorted_bam} > {output.flagstat}
        samtools stats -@ {threads} {input.sorted_bam} > {output.stats}
        """


##
## remove intermediate FASTQ and BAM files
##

rule remove_fastqs:
    output:
        remove_fastqs = "{outDir}/log/{sample}_removed_fastqs.txt"
    params:
        fastq1_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_1.fq.gz",
        fastq2_sortmerna = lambda wildcards: f"{outDir}/FASTQ_SORTMERNA/{wildcards.sample}_sortmerna_2.fq.gz",
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        unpairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_unpaired.fastq.gz",
        unpairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_unpaired.fastq.gz",
    conda:
        "../envs/bash.yaml"
    shell:
        """
        rm {params.fastq1_sortmerna} {params.fastq2_sortmerna}
        rm {params.pairedFile1} {params.pairedFile2} {params.unpairedFile1} {params.unpairedFile2}
        touch {output.remove_fastqs}
        """


rule remove_bams:
    input:
        f"resources/{date_str}_MULTIQC_end.txt"
    output:
        remove_bams = "{outDir}/log/{sample}_removed_bams.txt"
    params:
        aligned_bam = lambda wildcards: f"{outDir}/BAM/{wildcards.sample}_Aligned.out.bam",
        aligned_sorted_bam = lambda wildcards: f"{outDir}/BAM/{wildcards.sample}_Aligned.out.sorted.bam",
    conda:
        "../envs/bash.yaml"
    shell:
        """
        rm {params.aligned_bam}* {params.aligned_sorted_bam}*
        touch {output.remove_bams}
        """