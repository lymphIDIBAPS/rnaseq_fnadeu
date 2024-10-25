# RNA Seq Analysis Pipeline for F. Nadeu
![GitHub last commit](https://img.shields.io/github/last-commit/lymphIDIBAPS/rnaseq_fnadeu)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/lymphIDIBAPS/rnaseq_fnadeu)
----
This is a pipeline written in Python and bash, and with Snakemake as a workflow manager, that will output *kallisto* files (abundances of transcripts) from RNA-Seq data. These samples can be in .fastq or compressed format. 

The samples can be placed in any directory, but the path must be specified in the ***config/config.yaml*** file. 

# Table of Contents
1. [Installation](#installation)
2. [Linux](#linux)
3. [Windows](#windows)
4. [macOS](#macos)
5. [Snakemake Environment](#snakemake-environment)
6. [Clone the repository](#clone-the-repository)
7. [Snakemake Usage](#snakemake-usage)
8. [Run the pipeline in a HPC](#run-the-pipeline-in-a-hpc)
9. [Configuration of the pipeline](#configuration-of-the-pipeline)
10. [Cluster Configuration](#cluster-configuration)

## Installation


### Linux
<details><summary>Linux Installation</summary>

#### Install Miniconda 3
Open a Linux shell, then run these three commands to quickly and quietly download the latest 64-bit Linux miniconda 3 installer, rename it to a shorter file name, silently install, and then delete the installer.
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:
```bash
~/miniconda3/bin/conda init bash
```
```zsh
~/miniconda3/bin/conda init zsh
```
You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
</details>

---
### Windows
<details><summary>Windows Installation</summary>

Since Windows does not have access to the majority of packages we need in the pipeline, we need to install *Linux on Windows*, also known as *WSL*. In the ***Windows Power Shell***:
```powershell
wsl --install
# This command will install the Ubuntu distribution of Linux.
```
If you run into an issue during the installation process, please check [the installation section of the troubleshooting guide](https://learn.microsoft.com/en-us/windows/wsl/troubleshooting#installation-issues).

#### Set up your Linux user info

Once you have installed WSL, you will need to create a user account and password for your newly installed Linux distribution. See the Best practices for [setting up a WSL development environment guide to learn more](https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password).

#### Install Miniconda 3
Once you have a working shell in your *WSL*, run these three commands to quickly and quietly download the latest 64-bit Linux miniconda 3 installer, rename it to a shorter file name, silently install, and then delete the installer.
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:
```bash
~/miniconda3/bin/conda init bash
```
```zsh
~/miniconda3/bin/conda init zsh
```
You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
</details>

----
### macOS
<details><summary>macOS Installation</summary>

These four commands download the latest M1 version of the MacOS installer, rename it to a shorter file name, silently install, and then delete the installer:
```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:
```bash
~/miniconda3/bin/conda init bash
```
```zsh
~/miniconda3/bin/conda init zsh
```

You should see ```(base)``` in the command line prompt. This tells you that you’re in your base conda environment. To learn more about conda environments, see [Environments](https://docs.anaconda.com/working-with-conda/environments/).

Check for a good installation with:
```bash
conda --version
# conda 24.X.X

conda list
# outputs a list of packages installed in the current environment (base)
```
<!-- #### Install mamba
Now we will install ***mamba***, that is a fast-ligthweigth package manager aking to conda.

In case you are downloading from the Clinic network, you will have trouble with the SLL certificate. TO solve any problems, do the following:
1. You can usually get a copy by clicking on the padlock icon in your browser when visiting any https site, then click around to view certificate, and download in PEM format.
Then we will point conda to it in our system. 
```bash
conda config --set ssl_verify <pathToYourFile>.pem
```

2. Next we will install mamba in our base environment with:
```bash
conda install -n base -c conda-forge mamba
```
To check:
```bash
mamba --version
# mamba 1.X.X
# conda 24.X.X
``` -->
</details>

---
### Snakemake Environment
Now, with miniconda installed in our machine, we can create a new environment with snakemake installed:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

In case you are downloading from the Clinic network, you will have trouble with the SLL certificate. To solve any problems, do the following:
1. You can usually get a copy by clicking on the padlock icon in your browser when visiting any https site, then click around to view certificate, and download in PEM format.

2. Then we will point conda to it in our system. 
```bash
conda config --set ssl_verify <pathToYourFile>.pem
```

Once installed, we must activate and move into the snakemake environment with:
```bash
conda activate snakemake
snakemake --version
# 8.16.X
```
If at any time we want to exit the environment, we can with```conda deactivate```, and to get back in with ```conda activate snakemake```.
To see the packages we have currently installed in the environment, we can with```conda list```.
### Clone the repository

1. Above the list of files, click Code.

![Imgur](https://i.imgur.com/FiesQMl.png)

2. Copy the URL for the repository. To clone the repository using HTTPS, under "HTTPS", copy the link provided.
3. Open a Terminal.
4. Change the current working directory to the location where you want the cloned directory. For example, ```cd rna_seq_fnadeu```. Make sure that the directory exists before you move into it.
5. Type ```git clone git@github.com:lymphIDIBAPS/rnaseq_fnadeu.git```.
6. Press Enter to create your local clone.
```bash
git clone git@github.com:lymphIDIBAPS/rnaseq_fnadeu.git
> Cloning into `rna_seq_fnadeu`...
> remote: Counting objects: 10, done.
> remote: Compressing objects: 100% (8/8), done.
> remove: Total 10 (delta 1), reused 10 (delta 1)
> Unpacking objects: 100% (10/10), done.
```
## Snakemake Usage
When we have the cloned repository, we can proced and add our sample data to the FASTQ directory. This is not mandatory, as in ***config/config.yaml*** file we can edit and set any path to our sample data. 

In the same file we can edit the number of threads our computer has, so it will run adapted to the current resources we have available. 

The rulegraph for our pipeline at date 25/10 is the following:
![Imgur](https://imgur.com/fW5OXu6.png)

```bash
# For a test run of the pipeline
snakemake --use-conda -np

# For a real run of the pipeline
snakemake --use-conda
```

## Run the pipeline in a HPC
If we have many samples and our computer does not have enough computational power, we can run the pipeline in a cluster. This pipeline has been prepared to run in the [StarLife](https://www.bsc.es/supportkc/docs/StarLife/intro) cluster, in the [BSC](https://www.bsc.es/).

1. Make a new directory named ***/slgpfs/*** in your computer and mount it to the same directory in StarLife:
```bash
mkdir /home/user/slgpfs
sshfs -o allow_other your_bsc_user@sl1.bsc.es:/slgpfs/ /home/user/slgpfs/
```
This will allow you to see and work on the custer from your computer system directly.

2. On your computer, navigate to the directory: ***/home/user/slgpfs/projects/group_folder***

3. Download and extract the following file to the directory, in which we have a full conda environment ready to run snakemake:
[Snakemake Conda Environment](https://drive.google.com/file/d/1ZQB02jZhhr5G_DlbhDKFHiqZ7AVuCxqx/view?usp=sharing)

4. Clone this repository in the directory, following the steps from [Clone the repository](#clone-the-repository)

5. Now, connect to the cluster:
```bash
ssh your_username@sl1.bsc.es # or
ssh your_username@sl2.bsc.es
```
6. In the cluster, navigate to the cloned repository: ***/slgpfs/projects/group_folder/rna_seq_fnadeu***

7. Now, activate the snakemake_bsc environment:
```bash
source ../snakemake_bsc/bin/activate
```
In your terminal, you should now see something like: ```(snakemake_bsc) your_username@sllogin1```

8. At this point, in your local machine, you can move your samples to the directory ***/slgpfs/projects/group_folder/rna_seq_fnadeu/FASTQ***. 

9. Now, you can run the pipeline from the cluster with the command:
```bash
# For a test run of the pipeline
snakemake --profile config/slurm/ --use-envmodules -np

# For a real run of the pipeline
snakemake --profile config/slurm/ --use-envmodules
``` 
This command above will run the pipeline with the pipeline configuration from the file located in ***/slgpfs/projects/group_folder/rna_seq_fnadeu/config/config.yaml***. Be sure to check and modify the configuration file to alter the pipeline with your desired options.

The cluster configuration file is located in ***/home/oscar/rnaseq/config/slurm/config.yaml***. Below you have all the options available to customize your cluster run.

## Configuration of the pipeline
### General Configuration

1. Fastq directory paired unmerged: specify the directory containing the paired end unmerged FASTQ files

2. Technical duplicates: specify if we must work with tecnical duplicates or not

3. Sample extension: specify the suffix of the samples

4. Fastq directory single end: specify the directory containing the single end FASTQ files

5. Perform FASTQC: specify if you want to run FASTQC quality control on the samples

6. Perform multiQC: specify if you want to run multiQC quality control on the results from FASTQC, sortmerna and trimmomatic. Perform FASTQC must be set to yes.

### Trimmomatic
In this section of the configuration file we can adjust the options related with the trimmomatic tool.
The current defult values are the same as those used by Marta Sureda. 

#### ILUMINACLIP
This step is used to find and remove Illumina adapters.

1. Seed mismatches:  specifies the maximum mismatch count which will still allow a full
match to be performed.

2. Palindrome clip threshold: specifies how accurate the match between the two 'adapter
ligated' reads must be for PE palindrome read alignment.

3. Simple clip threshold: pecifies how accurate the match between any adapter sequence must be against a read.

#### SLIDINGWINDOW
Perform a sliding window trimming, cutting once the average quality within the window falls
below a threshold.
1. Window Size: specifies the number of bases to average across.

2. Required quality: specifies the average quality required.
#### LEADING
Remove low quality bases from the beginning. As long as a base has a value below this
threshold the base is removed and the next base will be investigated.

1. Leading quality: specifies the minimum quality required to keep a base.
#### TRAILING
Remove low quality bases from the end. As long as a base has a value below this threshold
the base is removed and the next base will be investigated.

1. Trailing quality: specifies the minimum quality required to keep a base.
#### MINLEN
This module removes reads that fall below the specified minimal length. If required, it should normally be after all other processing steps.

1. Minimal length: specifies the minimum length of reads to be kept.

### Kallisto Single End
For the single end mode you need to supply the ``--single`` flag as well as the ``-l`` and ``-s``options: 

``-l, --fragment-length=DOUBLE``: estimated average fragment length

``-s, --sd=DOUBLE``: estimated standard deviation of fragment length

## Cluster Configuration
Remember to check the files in ***/config/slurm/config.yaml*** for the cluster configuration. Review all the items and in case something is not clear you can check [in this website](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html#advanced-resource-specifications) what each term means in the configuration. 

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## Author
Developed by [@obaeza16](https://github.com/obaeza16), based on a pipeline written by Marta Sureda.

Mantained by [Lymphoid neoplasms program, IDIBAPS](https://www.clinicbarcelona.org/en/idibaps/programs/lymphoid-neoplasms-programme) for Ferran Nadeu. 

## License

[MIT](https://choosealicense.com/licenses/mit/)