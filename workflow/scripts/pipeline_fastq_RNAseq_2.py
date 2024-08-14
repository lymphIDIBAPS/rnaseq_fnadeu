# Libraries:
import subprocess 
import argparse
import time
import sys

# Input:
parser  = argparse.ArgumentParser(prog='pipeline_fastq_RNAseq_2', description='''RNA-seq pipeline (script 2)''')

parser.add_argument('-p', '--pathToReferenceDataAndPipelines',
                    dest = "pathToReferenceDataAndPipelines",
                    action = "store",
                    help = "Path to folder reference_data_and_pipelines")

parser.add_argument('-o', '--outDir',
                    dest = "outDir",
                    action = "store",
                    required=True,
                    help = "Working directory")

parser.add_argument('-s', '--sample',
                    dest = "sample",
                    action = "store",
                    help = "Sample name")

parser.add_argument('-F1', '--fastq1',
					dest = "fq1",
					action = "store",
					required=True,
					help = "Fastq file R1")

parser.add_argument('-F2', '--fastq2',
					dest = "fq2",
					action = "store",
					required=True,
					help = "Fastq file R2")

parser.add_argument('-TS', '--TranscriptionStrand',
					dest = "TranscriptionStrand",
					action = "store",
					required=True,
					help = "Transcription strand [first, second, unstranded]")

parser.add_argument('-tT', '--tTrimming',
					dest = "tTrimming",
					action = "store",
					default = "no",
					help = "Perform T-trimming")

parser.add_argument('-ad', '--adapters',
					dest = "adapters",
					action = "store",
					default = "illumina",
                    choices=['illumina', 'bioskryb'],
					help = "Adapters to be removed. Choices: illumina or bioskryb; default = illumina")

parser.add_argument('-r', '--removeFastqs',
                    dest = "removeFastqs",
                    action = "store",
                    default = "yes",
                    help = "Remove intermediate fastqs")

parser.add_argument('-rb', '--removeBams',
                    dest = "removeBams",
                    action = "store",
                    default = "yes",
                    choices=['yes', 'no'],
                    help = "Remove bam files after QC (yes or no, default=yes)")

parser.add_argument('-sortmernaDB', '--sortmernaDB',
                    dest = "sortmernaDB",
                    action = "store",
					default = "default",
                    choices=['fast', 'default', 'sensitive'],
                    help = "rRNA database for sortmerna")

parser.add_argument('-g', '--genes',
                    dest = "genes",
                    action = "store",
                    default = "cDNA",
                    choices=['cDNA', 'ncRNA', 'both'],
                    help = "Index file for kallisto including only cDNA, ncRNA, or both")

parser.add_argument('-qc', '--runQC',
                    dest = "runQC",
                    action = "store",
                    default = "yes",
                    choices=['yes', 'no'],
                    help = "Should BAM file be created and QC metrics analyzed? (default = yes)")

parser.add_argument('-@', '--cpus',
                    dest = "cpus",
                    action = "store",
                    help = "Number of cpus per task")

options = parser.parse_args()


# 0) define variables and files
sortmerna = options.pathToReferenceDataAndPipelines+"/programs/sortmerna-4.3.4/bin/sortmerna"
sortmerna_db = options.pathToReferenceDataAndPipelines+"/programs/sortmerna-4.3.4/db/smr_v4.3_"+options.sortmernaDB+"_db.fasta"
star_genome = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/STAR_100"
refFlat = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GFF3/Homo_sapiens.GRCh38.105.refFlat"
ribosomal_intervals = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GTF/Homo_sapiens.GRCh38.105.ribosomalIntervals"
rseqc_bed = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/ensembl_GFF3/Homo_sapiens.GRCh38.105.bed"
if options.genes == "cDNA": kallisto_ref = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.release-105.idx"
elif options.genes == "ncRNA": kallisto_ref = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.ncrna.release-105.idx"
else: kallisto_ref = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38.cdna.all.ncrna.release-105.idx"
outDir = options.outDir
sample = options.sample
fastq1 = options.fq1
fastq2 = options.fq2
TranscriptionStrand = options.TranscriptionStrand
cpus = options.cpus

start_time = time.time()
LOG_FILE_OUT = open(outDir+"/log/log_out_"+sample+".txt", 'w')
LOG_FILE_ERR = open(outDir+"/log/log_err_"+sample+".txt", 'w')


##
## 1) sortmerna
##
# create 'sample' dir
bashArguments = "mkdir "+outDir+"/FASTQ_SORTMERNA/"+sample
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

# run sortmerna
bashArguments = sortmerna+" --ref "+sortmerna_db+" --reads "+fastq1+" --reads "+fastq2+" --other --workdir "+outDir+"/FASTQ_SORTMERNA/"+sample+" --fastx --paired_in -threads "+cpus+" -out2 -v"
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
    sys.exit(1)

# move/adjust outputs
bashArguments = "rm -rf "+outDir+"/FASTQ_SORTMERNA/"+sample+"/idx/; rm -rf "+outDir+"/FASTQ_SORTMERNA/"+sample+"/kvdb/; rm -rf "+outDir+"/FASTQ_SORTMERNA/"+sample+"/readb/"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

bashArguments = "mv "+outDir+"/FASTQ_SORTMERNA/"+sample+"/out/other_fwd.fq.gz "+outDir+"/FASTQ_SORTMERNA/"+sample+"_sortmerna_1.fq.gz"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

bashArguments = "mv "+outDir+"/FASTQ_SORTMERNA/"+sample+"/out/other_rev.fq.gz "+outDir+"/FASTQ_SORTMERNA/"+sample+"_sortmerna_2.fq.gz"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

bashArguments = "mv "+outDir+"/FASTQ_SORTMERNA/"+sample+"/out/aligned.log "+outDir+"/FASTQ_SORTMERNA/"+sample+"_aligned.log"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

logInput = open(outDir+"/FASTQ_SORTMERNA/"+sample+"_aligned.log", "r")
logSortmerna = open(outDir+"/MULTIQC/files/"+sample+"_aligned.log", "w")
for lLine in logInput:
	logSortmerna.write(lLine.replace(fastq1.split("/")[-1], sample+".fq.gz").replace(fastq2.split("/")[-1], sample+".fq.gz"))
logInput.close()
logSortmerna.close()

bashArguments = "rm -rf "+outDir+"/FASTQ_SORTMERNA/"+sample+"/out/; rm -rf "+outDir+"/FASTQ_SORTMERNA/"+sample
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)



##
## 2) trimming
##
# prepare paths to inputs and outputs
fastq1_sortmerna = outDir+"/FASTQ_SORTMERNA/"+sample+"_sortmerna_1.fq.gz"
fastq2_sortmerna = outDir+"/FASTQ_SORTMERNA/"+sample+"_sortmerna_2.fq.gz"
pairedFile1 = outDir+"/FASTQ_TRIMMED/"+sample+"_sortmerna_1_paired.fastq.gz"
unpairedFile1 = outDir+"/FASTQ_TRIMMED/"+sample+"_sortmerna_1_unpaired.fastq.gz"
pairedFile2 = outDir+"/FASTQ_TRIMMED/"+sample+"_sortmerna_2_paired.fastq.gz"
unpairedFile2 = outDir+"/FASTQ_TRIMMED/"+sample+"_sortmerna_2_unpaired.fastq.gz"
if options.adapters == "illumina":
	adaptersSeq = "/slgpfs/projects/cli79/reference_data_and_pipelines/programs/TRIMMOMATIC/adapters/TruSeq3-PE-2.fa"
elif options.adapters == "bioskryb":
	adaptersSeq = "/slgpfs/projects/cli79/reference_data_and_pipelines/programs/TRIMMOMATIC/adapters/Bioskryb_ResolveOME.fa"

# run trimmomatic
if options.tTrimming == "no":
	bashArguments = "java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar PE -threads "+cpus+" -phred33 "+fastq1_sortmerna+" "+fastq2_sortmerna+" "+pairedFile1+" "+unpairedFile1+" "+pairedFile2+" "+unpairedFile2+" ILLUMINACLIP:"+adaptersSeq+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> "+outDir+"/FASTQ_TRIMMED/"+sample+"-trimmomatic.log"
else:
	bashArguments = "java -jar /apps/TRIMMOMATIC/0.40/trimmomatic-0.40-rc1.jar PE -threads "+cpus+" -phred33 "+fastq1_sortmerna+" "+fastq2_sortmerna+" "+pairedFile1+" "+unpairedFile1+" "+pairedFile2+" "+unpairedFile2+" ILLUMINACLIP:"+adaptersSeq+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:1 2> "+outDir+"/FASTQ_TRIMMED/"+sample+"-trimmomatic.log"
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
    sys.exit(1)

# move/adjust outputs
logInput = open(outDir+"/FASTQ_TRIMMED/"+sample+"-trimmomatic.log", "r")
logTrimmomatic = open(outDir+"/MULTIQC/files/"+sample+"-trimmomatic.log", "w")
for lLine in logInput:
	logTrimmomatic.write(lLine.replace("_sortmerna_1", "").replace("_sortmerna_2", ""))
logInput.close()
logTrimmomatic.close()



##
## 3) fastqc
##
bashArguments = "fastqc -o "+outDir+"/MULTIQC_FASTQC/files/ -t 6 "+fastq1+" "+fastq2+" "+fastq1_sortmerna+" "+fastq2_sortmerna+" "+pairedFile1+" "+pairedFile2
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
    sys.exit(1)



##
## 4) kallisto
##
# create 'sample' dir
bashArguments = "mkdir "+outDir+"/KALLISTO/"+sample
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

# open file to store stdout kallisto to use in multiqc
LOG_FILE_ERR_kallisto = open(outDir+"/KALLISTO/"+sample+"/log_err.txt", 'w')

# run kallisto
strand = " --rf-stranded" if TranscriptionStrand == "second" else " --fr-stranded" if TranscriptionStrand == "first" else ""
bashArguments = "kallisto quant -t "+cpus+" -i "+kallisto_ref+" -o "+outDir+"/KALLISTO/"+sample+" "+pairedFile1+" "+pairedFile2+strand
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR_kallisto)
if e != 0:
    LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
    LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
    sys.exit(1)

# close file to store stdout
LOG_FILE_ERR_kallisto.close()

# move/adjust outputs
bashArguments = "mv "+outDir+"/KALLISTO/"+sample+"/abundance.h5 "+outDir+"/KALLISTO/"+sample+"_abundance.h5"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

bashArguments = "mv "+outDir+"/KALLISTO/"+sample+"/abundance.tsv "+outDir+"/KALLISTO/"+sample+"_abundance.tsv"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

bashArguments = "mv "+outDir+"/KALLISTO/"+sample+"/run_info.json "+outDir+"/KALLISTO/"+sample+"_run_info.json"
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

LOG_FILE_ERR_kallisto = open(outDir+"/KALLISTO/"+sample+"/log_err.txt", 'r')
logKallistoMultiqc = open(outDir+"/MULTIQC/files/"+sample+"_log_err.txt", 'w')
for lLine in LOG_FILE_ERR_kallisto:
	LOG_FILE_ERR.write(lLine)
	logKallistoMultiqc.write(lLine.replace("_sortmerna_1_paired", "").replace("_sortmerna_2_paired", ""))
logKallistoMultiqc.close()

bashArguments = "rm -rf "+outDir+"/KALLISTO/"+sample
subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)


###
### BAM + QC if specified
###
if options.runQC == "yes":
    ##
    ## 5) align
    ##
    # run STAR
    bashArguments = "STAR --runThreadN "+cpus+" --genomeDir "+star_genome+" --readFilesIn "+pairedFile1+" "+pairedFile2+" --readFilesCommand zcat --outFileNamePrefix "+outDir+"/BAM/"+sample+"_ --outSAMtype BAM Unsorted --twopassMode Basic --outBAMcompression 10" # not used now compared to Romina's code: --outSAMstrandField intronMotif   
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # sort using samtools
    bashArguments = "samtools sort -@ "+cpus+" -m 3G -o "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam "+outDir+"/BAM/"+sample+"_Aligned.out.bam"
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # index sorted bam using samtools
    bashArguments = "samtools index -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam"
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # move/adjust outputs
    bashArguments = "rm -rf "+outDir+"/BAM/"+sample+"__STARgenome/; rm -rf "+outDir+"/BAM/"+sample+"__STARpass1/"
    subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

    bashArguments = "cp "+outDir+"/BAM/"+sample+"_Log.final.out "+outDir+"/MULTIQC/files/"+sample+"_Log.final.out"
    subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

    # md5sum
    subprocess.call("md5sum "+outDir+"/BAM/"+sample+"_Aligned.out.bam > "+outDir+"/BAM/"+sample+"_Aligned.out.bam.md5", shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    subprocess.call("md5sum "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam.md5", shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    subprocess.call("md5sum "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam.bai > "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam.bai.md5", shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)



    ##
    ## 6) QC
    ##
    # CollectRnaSeqMetrics
    strand = "SECOND_READ_TRANSCRIPTION_STRAND" if TranscriptionStrand == "second" else "FIRST_READ_TRANSCRIPTION_STRAND" if TranscriptionStrand == "first" else "NONE"
    bashArguments = "java -Xmx"+str(int(cpus)*4)+"G -jar /apps/PICARD/2.24.0/bin/picard.jar CollectRnaSeqMetrics -I "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam "+" -O "+outDir+"/MULTIQC/files/"+sample+".CollectRnaSeqMetrics -REF_FLAT "+refFlat+" -STRAND "+strand+" -RIBOSOMAL_INTERVALS "+ribosomal_intervals
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # RSeQC read_distribution.py
    bashArguments = "read_distribution.py -i "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam -r "+rseqc_bed+" > "+outDir+"/MULTIQC/files/"+sample
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # # RSeQC geneBody_coverage.py ==> it takes a lot of time to run compared to the entire pipeline (60% of execution time) and we get a similar output with CollectRnaSeqMetrics
    # bashArguments = "geneBody_coverage.py -i "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam -r "+rseqc_bed+" -o "+outDir+"/MULTIQC/files/"+sample
    # e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    # if e != 0:
    #	LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
    #	LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
    #	sys.exit(1)

    # samtools idxstats
    bashArguments = "samtools idxstats -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".idxstats"
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # samtools flagstat
    bashArguments = "samtools flagstat -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".flagstat"
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)

    # samtools stats
    bashArguments = "samtools stats -@ "+cpus+" "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam > "+outDir+"/MULTIQC/files/"+sample+".stats"
    e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    if e != 0:
        LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
        LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
        sys.exit(1)



##
## 7) remove intermediate FASTQ and BAM files
##
if options.removeFastqs == "yes":
	# sortmerna
	bashArguments = "rm "+fastq1_sortmerna+" "+fastq2_sortmerna
	subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)

	# trimmomatic
	bashArguments = "rm "+pairedFile1+" "+unpairedFile1+" "+pairedFile2+" "+unpairedFile2
	subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
    
if options.runQC == "yes" and options.removeBams == "yes" :
	bashArguments = "rm "+outDir+"/BAM/"+sample+"_Aligned.out.bam* "+outDir+"/BAM/"+sample+"_Aligned.out.sorted.bam*"
	subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)



##
## DONE!
##
LOG_FILE_OUT.write("DONE. EXECUTION TIME IN MINUTES: %s\n" %(round((time.time() - start_time)/60, 4)))
LOG_FILE_OUT.close()
LOG_FILE_ERR.close()
