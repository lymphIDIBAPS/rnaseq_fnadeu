# Libraries:
import argparse
import os
import subprocess
import time
import glob
import sys

# Input:
parser  = argparse.ArgumentParser(prog='pipeline_fastq_RNAseq_1', description='''RNA-seq pipeline (script 1)''')

parser.add_argument('-p', '--pathToReferenceDataAndPipelines',
                    dest = "pathToReferenceDataAndPipelines",
                    action = "store",
                    default = "/slgpfs/projects/cli79/reference_data_and_pipelines",
                    help = "Path to folder reference_data_and_pipelines")

parser.add_argument('-i', '--infoRun',
                    dest = "infoRun",
                    required=True,
                    action = "store",
                    help = "Path to RunInfo.txt containing Sample\tSeq\tFastq1\Fastq2\tTranscriptionStrand\t... plus additional metadata")

parser.add_argument('-w', '--workDir',
                    dest = "workDir",
                    required=True,
                    action = "store",
                    help = "Working directory")

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

parser.add_argument('-t', '--time',
                    dest = "time",
                    action = "store",
                    default = "72",
                    help = "Maximum execution time (in hours)")

parser.add_argument('-n', '--analysisName',
                    dest = "analysisName",
                    action = "store",
                    default = "",
                    help = "Analysis name")

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
                    help = "Remove bam files after QC (yes or no, default = yes)")

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
                    help = "Index file for kallisto including only cDNA, ncRNA, or both (default = cDNA)")

parser.add_argument('-qc', '--runQC',
                    dest = "runQC",
                    action = "store",
                    default = "yes",
                    choices=['yes', 'no'],
                    help = "Should BAM file be created and QC metrics analyzed? (default = yes)")

parser.add_argument('-@', '--cpus',
                    dest = "cpus",
                    action = "store",
                    default = "20",
                    help = "Number of cpus per job (default = 20)")

options = parser.parse_args()

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

##
## Prepare and submit job script for each sample
##
nTasks = 0
jobIds = list()
I = open(options.infoRun, "r")
for iLine in I:

    if iLine.startswith("Sample\t") or iLine.startswith("#"): continue
    
    iList = iLine.rstrip("\n").split("\t")
    sample = iList[0]
    fq1 = iList[2]
    fq2 = iList[3]
    strand = iList[4]

    # Prepare command
    comm = "python "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_2.py -p "+options.pathToReferenceDataAndPipelines+" -o "+outDir+" -s "+sample+" -F1 "+fq1+" -F2 "+fq2+" -TS "+strand+" -tT "+options.tTrimming+" -ad "+options.adapters+" -r "+options.removeFastqs+" -rb "+options.removeBams+" -sortmernaDB "+options.sortmernaDB+" -g "+options.genes+" -qc "+options.runQC+" -@ "+options.cpus
    
    # Write job script
    O = open(outDir+"/job_script_"+sample+".sh", "w")
    O.write("#!/bin/bash\n")
    O.write("#SBATCH --job-name=RNAseq_"+sample+"\n")
    O.write("#SBATCH -D '.'\n")
    O.write("#SBATCH --output=%j.out\n")
    O.write("#SBATCH --error=%j.err\n")
    O.write("#SBATCH --ntasks=1\n") 
    O.write("#SBATCH --time="+options.time+":00:00\n")
    O.write("#SBATCH --cpus-per-task="+options.cpus+"\n")
    if int(options.time) <= 2: O.write("#SBATCH --qos=debug\n")
    O.write("\nmodule load python/3.6.5 java/10.0.1 fastqc/0.11.9 gcc/9.2.0 zlib/1.2.11 hdf5/1.10.2 szip/2.1.1 kallisto/0.46.1 star/2.7.8a samtools/1.9 picard/2.24.0 R/3.5.0\n\n") 
    O.write(comm+"\n")
    O.close()

    # Run job and store jobId
    os.chdir(outDir)
    process = subprocess.Popen("sbatch job_script_"+sample+".sh", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    print("Submitted batch job "+output.decode("utf-8").rstrip("\n").split(" ")[-1]+" ---> sample: "+sample)
    jobIds.append(output.decode("utf-8").rstrip("\n").split(" ")[-1])
    os.chdir(currentDir)

    nTasks += 1

I.close()


##
## Summary QC and plots once all previous jobs are done
##
# multiqc
comm1 = "multiqc -f -i "+date_str+"_pipeline_fastq_RNAseq_FASTQC"+aName+" -b 'Multiqc report for RNAseq pipeline (FASTQC)' -n "+date_str+"_pipeline_fastq_RNAseq_FASTQC"+aName+" -o "+outDir+"/MULTIQC_FASTQC "+outDir+"/MULTIQC_FASTQC/files"    
comm2 = "multiqc -f -i "+date_str+"_pipeline_fastq_RNAseq_BAM"+aName+" -b 'Multiqc report for RNAseq pipeline' -n "+date_str+"_pipeline_fastq_RNAseq_BAM"+aName+" -o "+outDir+"/MULTIQC "+outDir+"/MULTIQC/files"

# plots
ensemblTable = options.pathToReferenceDataAndPipelines+"/data/genome_GRCh38.p13_GCA_000001405.28/kallisto/Homo_sapiens.GRCh38_genes_transcripts_release-105.txt"
comm3 = "Rscript --vanilla "+options.pathToReferenceDataAndPipelines+"/RNAseq/pipeline_fastq_RNAseq_3.R "+ensemblTable+" "+options.infoRun+" "+outDir+"/KALLISTO "+outDir+"/MULTIQC/"+date_str+"_pipeline_fastq_RNAseq_PCAs.pdf"

# Prepare job script
O = open(outDir+"/job_script_summary.sh", "w")
O.write("#!/bin/bash\n")
O.write("#SBATCH --job-name=RNAseq_summary\n")
O.write("#SBATCH -D '.'\n")
O.write("#SBATCH --output=%j.out\n")
O.write("#SBATCH --error=%j.err\n")
O.write("#SBATCH --ntasks=1\n") 
O.write("#SBATCH --time="+options.time+":00:00\n")
O.write("#SBATCH --cpus-per-task=2\n")
if int(options.time) <= 2: O.write("#SBATCH --qos=debug\n")
O.write("\nmodule load python/3.6.5 java/10.0.1 fastqc/0.11.9 gcc/9.2.0 zlib/1.2.11 hdf5/1.10.2 szip/2.1.1 kallisto/0.46.1 star/2.7.8a samtools/1.9 picard/2.24.0 R/3.5.0\n\n") 
O.write("export LC_ALL=C.UTF-8\n")
O.write("export LANG=C.UTF-8\n\n")
O.write(comm1+"\n")
O.write(comm2+"\n")
O.write(comm3+"\n")
O.close()

# Run job
os.chdir(outDir)
dependencies = ":".join(jobIds)
print("Submitting summary script...")
subprocess.call("sbatch --dependency=afterany:"+dependencies+" job_script_summary.sh", shell=True)
nTasks += 1


# Print number of jobs submitted
print("Number of jobs submitted: "+str(nTasks))