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
THREADS = int(config["cpus"])

if config["genes"] == "cDNA": 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.release-105.idx"
elif config["genes"] == "ncRNA": 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.ncrna.release-105.idx"
else: 
    kallisto_ref = config["pathToReferenceDataAndPipelines"]+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.ncrna.release-105.idx"

# Extract sample names and paths to reads
samples = pep.sample_table

# Functions to extract flag values

def check_trimming(tTrimming):
    return "HEADCROP:1" if tTrimming == "yes" else ""

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

# Module loads
"python/3.6.5" "java/10.0.1" "fastqc/0.11.9" 
"gcc/9.2.0" "zlib/1.2.11" "hdf5/1.10.2" "szip/2.1.1" "kallisto/0.46.1" 
"star/2.7.8a" "samtools/1.9" "picard/2.24.0" "R/3.5.0"

# Pipeline Rules

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
    threads:
        THREADS
    envmodules:
        "python/3.6.5"
        "sortmerna/4.3.6"
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

rule trimmomatic:
    input:
        config_file = "config/project_config.yaml"
    output:
        log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}_trimmomatic.log"
    params:
        fastq1=lambda wildcards: samples.loc[wildcards.sample]["forward"],
        fastq2=lambda wildcards: samples.loc[wildcards.sample]["reverse"],
        paired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        unpaired1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_unpaired.fastq.gz",
        paired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        unpaired2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_unpaired.fastq.gz",
        adapters = "resources/TruSeq3-PE.fa",  # Default adapter file
        result = check_trimming(config["tTrimming"])   # Default trimming option
    threads:
        THREADS
    envmodules:
        "java/10.0.1"
    shell:
        """
        java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar PE -threads {threads} -phred33 {params.fastq1} {params.fastq2} \
        {params.paired1} {params.unpaired1} {params.paired2} {params.unpaired2} \
        ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 {params.result} 2> {output.log}
        """

# java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar

rule multiqc_trimmomatic_log:
    input:
        aligned_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}_trimmomatic.log",
    output:
        multiqc_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/MULTIQC/files/{sample}_trimmomatic.log"
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
        aligned_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/FASTQ_TRIMMED/{sample}_trimmomatic.log",
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
    threads:
        THREADS
    envmodules:
        "python/3.6.5"
        "fastqc/0.11.9" 
    shell:
        """
        fastqc -o {params.outDir}/MULTIQC_FASTQC/files/ -t {threads} {params.fastq1} {params.fastq2} \
        {params.fastq1_sortmerna} {params.fastq2_sortmerna} {params.pairedFile1} {params.pairedFile2}
        """


## kallisto_index: builds an index
## from a FASTA formatted file of target sequences. Compute intensive rule

rule kallisto_index:
    input:
        index_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    output:
        index_out_path = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"
    threads:
        THREADS
    envmodules:
        "kallisto/0.46.1"
    shell:
        """
        kallisto index -i {output.index_out_path} --threads={threads} {input.index_path}
        """


rule kallisto:
    # This rule does not work with the paired-end file samples
    input:
        "resources/start.txt",
        "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx"
    output: 
        "{outDir}/KALLISTO/{sample}_abundance.h5",
        log = "{outDir}/KALLISTO/{sample}_kallisto.log"
    params:
        fastq1 = lambda wildcards: samples.loc[wildcards.sample]["forward"],
        transcription_strand = lambda wildcards: samples.loc[wildcards.sample]["TranscriptionStrand"],
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        kallisto_ref = "resources/kallisto/Homo_sapiens.GRCh38.cdna.all.release-100.idx",
    threads:
        THREADS
    run:
        # Get the strand flag based on the transcription strand
        # From the kallisto quant I deleted the {strand_flag} because it did not find any pseudoalignment with it
        strand_flag = get_strand_flag(params.transcription_strand)
        
        shell ("""
        mkdir -p {params.outDir}/KALLISTO/{wildcards.sample}
        touch {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5
        kallisto quant -t {threads} -i {params.kallisto_ref} -o {params.outDir}/KALLISTO/{wildcards.sample} {params.pairedFile1} {params.pairedFile2} 2> {output.log}
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.h5 {params.outDir}/KALLISTO/{wildcards.sample}_abundance.h5
        mv {params.outDir}/KALLISTO/{wildcards.sample}/abundance.tsv {params.outDir}/KALLISTO/{wildcards.sample}_abundance.tsv
        mv {params.outDir}/KALLISTO/{wildcards.sample}/run_info.json {params.outDir}/KALLISTO/{wildcards.sample}_run_info.json
        rm -rf {params.outDir}/KALLISTO/{wildcards.sample}
        """)

rule multiqc_kallisto_log:
    input:
        kallisto_log = "/home/oscar/RNAseq_ferran/20240902_162657_pipeline_fastq_RNAseq_TEST/KALLISTO/{sample}_kallisto.log"
    output:
        multiqc_log = "{outDir}/MULTIQC/files/{sample}_kallisto.log"
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


rule star_genome:
    # I touch the file start_align.txt
    input: 
        "resources/start_align.txt"
    output: 
        "resources/STAR/transcriptInfo.tab"
    params:
        fastaFiles = "resources/STAR/GRCh38.primary_assembly.genome.fa.chr19",
        gtfFile = "resources/STAR/gencode.v29.primary_assembly.annotation.chr19.gtf",
        genomeDir = "resources/STAR/",
    threads:
        THREADS
    envmodules:
        "python/3.6.5"
        "star/2.7.8a"
    shell:
        """
        mkdir -p resources/STAR
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genomeDir} --genomeFastaFiles {params.fastaFiles} --sjdbGTFfile {params.gtfFile} --genomeSAindexNbases 11
        """

rule star_map:
    input:
        "resources/STAR/transcriptInfo.tab"
    output: 
        "{outDir}/BAM/{sample}_Aligned.out.bam"
    params:
        pairedFile1 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_1_paired.fastq.gz",
        pairedFile2 = lambda wildcards: f"{outDir}/FASTQ_TRIMMED/{wildcards.sample}_sortmerna_2_paired.fastq.gz",
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        genomeDir = "resources/STAR",
    threads:
        THREADS
    envmodules:
        "python/3.6.5"
        "star/2.7.8a"
    # log:
    #     "{outDir}/BAM/{sample}_Log.out"
    shell:
        """
        mkdir -p resources/STAR
        STAR --runThreadN {threads} --genomeDir {params.genomeDir} --readFilesIn {params.pairedFile1} {params.pairedFile2} \
        --readFilesCommand zcat --outFileNamePrefix {params.outDir}/BAM/{wildcards.sample}_ --outSAMtype BAM Unsorted --twopassMode Basic --outBAMcompression 10
        """
# bashArguments = "STAR --runThreadN "+cpus+" --genomeDir "+star_genome+" --readFilesIn "+pairedFile1+" "+pairedFile2+" --readFilesCommand zcat 
# --outFileNamePrefix "+outDir+"/BAM/"+sample+"_ --outSAMtype BAM Unsorted --twopassMode Basic --outBAMcompression 10" # not used now compared to Romina's code: --outSAMstrandField intronMotif   


rule samtools:
    input:
        unsorted_bam = "{outDir}/BAM/{sample}_Aligned.out.bam",
    output:
        sorted_bam = "{outDir}/BAM/{sample}_Aligned.out.sorted.bam",
    params:
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        star_genome = "{outDir}/BAM/{sample}__STARgenome/",
        star_pass1 = "{outDir}/BAM/{sample}__STARpass1/",
        log_final_out = "{outDir}/BAM/{sample}_Log.final.out",
        log_final_multiqc = "{outDir}/MULTIQC/files/{sample}_Log.final.out",
    threads:
        THREADS
    envmodules:
        "python/3.6.5"
        "samtools/1.9"
    shell:
        """
        samtools sort -@ {threads} -m 3G -o {output.sorted_bam} {input.unsorted_bam}
        samtools index -@ {threads} {output.sorted_bam}
        rm -rf {params.star_genome} {params.star_pass1}
        cp {params.log_final_out} {params.log_final_multiqc}
        """

# bashArguments = "samtools sort -@ "+cpus+" -m 3G -o "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam "+outDir+"/BAM/"+sample+"_Aligned.out.bam"
# bashArguments = "samtools index -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam"

# bashArguments = "rm -rf "+outDir+"/BAM/"+sample+"__STARgenome/; rm -rf "+outDir+"/BAM/"+sample+"__STARpass1/"
# bashArguments = "cp "+outDir+"/BAM/"+sample+"_Log.final.out "+outDir+"/MULTIQC/files/"+sample+"_Log.final.out"


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
        transcription_strand = lambda wildcards: samples.loc[wildcards.sample]["TranscriptionStrand"],
        outDir = config["workDir"] + "/" + "20240902_162657_pipeline_fastq_RNAseq" + aName,
        refFlat = "resources/refFlat.txt",
        ribosomal_intervals = "resources/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_v0_annotation_gencode_v34_rRNA.interval_list"
    run:
        # Get the strand flag based on the transcription strand
        # From the kallisto quant I deleted the {strand_flag} because it did not find any pseudoalignment with it
        strand_flag = get_transcription_strand(params.transcription_strand)

        shell ("""
        module load picard/2.24.0
        picard CollectRnaSeqMetrics -I {params.outDir}/BAM/{wildcards.sample}_Aligned.out.sorted.bam \
        -O {params.outDir}/MULTIQC/files/{wildcards.sample}.CollectRnaSeqMetrics -REF_FLAT {params.refFlat} -STRAND {strand_flag} -RIBOSOMAL_INTERVALS {params.ribosomal_intervals}
        """)

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
        "python/3.6.5"
        "samtools/1.9"
    shell:
        """
        samtools idxstats -@ {threads} {input.sorted_bam} > {output.idxstats}
        samtools flagstat -@ {threads} {input.sorted_bam} > {output.flagstat}
        samtools stats -@ {threads} {input.sorted_bam} > {output.stats}
        """

# bashArguments = "samtools idxstats -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".idxstats"
# bashArguments = "samtools flagstat -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".flagstat"
# bashArguments = "samtools stats -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".stats"


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
    shell:
        """
        echo {params.fastq1_sortmerna} {params.fastq2_sortmerna}
        echo {params.pairedFile1} {params.pairedFile2} {params.unpairedFile1} {params.unpairedFile2}
        touch {output.remove_fastqs}
        """


rule remove_bams:
    output:
        remove_bams = "{outDir}/log/{sample}_removed_bams.txt"
    params:
        aligned_bam = lambda wildcards: f"{outDir}/BAM/{wildcards.sample}_Aligned.out.bam",
        aligned_sorted_bam = lambda wildcards: f"{outDir}/BAM/{wildcards.sample}_Aligned.out.sorted.bam",
    threads:
        THREADS
    shell:
        """
        echo {params.aligned_bam} {params.aligned_sorted_bam}
        touch {output.remove_bams}
        echo {threads}
        """
# bashArguments = "rm "+outDir+"/BAM/"+sample+"_Aligned.out.bam* "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam*"

### DONE WITH RNAseq_2.py