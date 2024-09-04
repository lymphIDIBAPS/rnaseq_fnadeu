import os
import yaml
# import subprocess
import time
import peppy

with open("config/config.yaml", "r") as file:
    config = yaml.safe_load(file)

date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")

# Get current directory
currentDir = os.getcwd()

# Create output directory
aName = "_" + config["analysisName"] if config["analysisName"] != "" else ""
outDir = config["workDir"] + "/" + date_str + "_pipeline_fastq_RNAseq" + aName
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



# Load the PEP file
pep = peppy.Project("config/project_config.yaml")

# Extract sample names and paths to reads
samples = pep.sample_table

# Extract samples paths to reads
samples2 = pep.sample_table["sample_name"]

# Function to look up forward and reverse reads

### GPT

# Function to retrieve the correct file based on sample and fastq number
def get_fastq_file(sample, fastq):
    # Locate the row based on the sample name
    sample_row = samples.loc[sample]
    
    # Get the index of the desired fastq number
    try:
        fastq_index = sample_row['fastq'].index(fastq)
    except ValueError:
        raise ValueError(f"Fastq number {fastq} not found for sample {sample}.")
    
    # Return the file corresponding to the fastq index
    return sample_row['file'][fastq_index]

# Example usage:
file_name = get_fastq_file("human_3", "1")
print(file_name)


forw = get_fastq_file("paired-end", "1")
reve = get_fastq_file("paired-end", "2")


# KALLISTO PREPARATIONS

def get_strand_flag(transcription_strand):
    if transcription_strand == "first":
        return " --fr-stranded"
    elif transcription_strand == "second":
        return " --rf-stranded"
    else:
        return ""


strand_set = get_strand_flag(samples.loc["paired-end"]["TranscriptionStrand"])


threads = int(config["cpus"])
type(threads)