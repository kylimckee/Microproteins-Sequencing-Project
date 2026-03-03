# Microproteins-Sequencing-Project

Pipeline used to analyze sarcoma bulk RNA-sequencing data to identify novel open reading frames (ORFs) with the potential to encode microproteins.

---

## Software Setup

You will need to install Miniforge into our data directory using the method described in the [Conda on Biowulf](https://hpc.nih.gov/docs/diy_installation/conda.html) documentation.


### Change Working Directory 

```bash
cd /data/mckeeka
pwd
```

### Install MiniConda 

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
#Type "yes" to accept the license. When asked for the installation location, change the directory to "/data/mckeeka/miniconda3". Type "yes" when asked to initialize MiniConda.
```

### Activate Changes and Source bashrc File

```bash
source ~/.bashrc
conda --version
```

## Upload Data and Reference Into Directory

### Make Working Directory

```bash
mkdir bulkRNA_sarcoma
cd /data/mckeeka/bulkRNA_sarcoma
```

### Download and Reformat FASTQ Files

Sarcoma FASTQ files were transferred to a folder called "MCI_fastq_117_STS_FASTQ" in the bulkRNA_sarcoma directory using Globus.
The following code changes the file names to fit the pipeline format (<sample>.fastq.<read>.gz).

```bash
cd /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ
for f in *R1.fastq.gz; do
  sample=$(echo "$f" | sed -E 's/(.*)\.R1\.fastq\.gz/\1/')
  mv "$f" "${sample}.fastq.1.gz"
done

for f in *R2.fastq.gz; do
  sample=$(echo "$f" | sed -E 's/(.*)\.R2\.fastq\.gz/\1/')
  mv "$f" "${sample}.fastq.2.gz"
done
```

### Download Human Reference Transcriptome

The human reference transcriptome (Ensembl release 115, GRCh38) was downloaded into a reference directory within bulkRNA_sarcoma.

```bash
cd /data/mckeeka/bulkRNA_sarcoma
mkdir reference
cd /data/mckeeka/bulkRNA_sarcoma/reference
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

### Visualize Working Directory 

Visualize the working directory to ensure proper setup.

```bash
cd ..
ls -R
```

## Create QC Pipeline Working Directory

The QC pipeline requires a working directory where the FASTQ files can be accessed. You can symlink these files instead of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir rawQC
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/rawQC
mkdir fastqc
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA

ln -s /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ/
```

## Generate QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after sequencing.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
conda create -n rawQC -c bioconda snakemake fastqc multiqc -y
conda activate rawQC
```

### Create Snakemake Raw QC Configuration File

```bash
nano rawQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("MCI_fastq_117_STS_FASTQ/{sample}.fastq.1.gz").sample

rule all:
    input:
        expand("rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip", sample=SAMPLES, read=[1,2]),
        "rawQC/rawQC_multiqc_report.html"

rule fastqc:
  input:
    "MCI_fastq_117_STS_FASTQ/{sample}.fastq.{read}.gz"
  output:
    html="rawQC/fastqc/{sample}.fastq.{read}_fastqc.html",
    zip="rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip"
  threads: 4
  shell:
    """
    fastqc -t {threads} -o rawQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip",
            sample=SAMPLES,
            read=[1,2])
  output:
    "rawQC/rawQC_multiqc_report.html"
  shell:
    """
    multiqc rawQC/fastqc -o run_bulkRNA/rawQC
    """
```

### Run Raw QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
sbatch --cpus-per-task=4 --mem=16G --time=06-00:00:00 \--wrap "snakemake -s rawQC_pipeline.smk -j 4"
```

## Create CutAdapt Pipeline Working Directory

The CutAdapt pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir CutAdapt
mkdir trimmed_FASTQ
mkdir logs
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/logs
mkdir logs_CutAdapt
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
```

## Generate CutAdapt Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.

### Install CutAdapt Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
conda create -n CutAdapt -c bioconda snakemake cutadapt -y
conda activate CutAdapt
```

### Create Snakemake CutAdapt Configuration File

```bash
nano CutAdapt_pipeline.smk

# Add the following code to the configuration file:

adapter = "AGATCGGAAGAG"     #Illumina Universal Adapter
minimum_length = 15          #Decreased from the recommended 20 since I am interested in smORFs
quality_trimming = "20,20"   #Recommended value
overlap = 5                  #Recommended value
threads = 4 

SAMPLES = glob_wildcards("MCI_fastq_117_STS_FASTQ/{sample}.fastq.1.gz").sample

rule all:
    input:
        expand("trimmed_FASTQ/{sample}.fastq.1.trimmed.gz", sample=SAMPLES),
        expand("trimmed_FASTQ/{sample}.fastq.2.trimmed.gz", sample=SAMPLES)

rule cutadapt_pe:
  input:
    r1 = "MCI_fastq_117_STS_FASTQ/{sample}.fastq.1.gz",
    r2 = "MCI_fastq_117_STS_FASTQ/{sample}.fastq.2.gz"
  output:
    r1 = "trimmed_FASTQ/{sample}.fastq.1.trimmed.gz",
    r2 = "trimmed_FASTQ/{sample}.fastq.2.trimmed.gz"
  log:
    "logs/logs_CutAdapt/{sample}.CutAdapt.log"
  shell:
    """
    cutadapt \
        -a {adapter} \
        -A {adapter} \
        -m {minimum_length} \
        -q {quality_trimming} \
        -O {overlap} \
        --pair-filter=any \
        --cores {threads} \
        -o {output.r1} \
        -p {output.r2} \
        {input.r1} {input.r2} > {log} 2>&1
    """

```

### Run CutAdapt Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
sbatch --cpus-per-task=4 --mem=16G --time=06-00:00:00 \--wrap "snakemake -s CutAdapt_pipeline.smk -j 4"
```

## Create Indexing Pipeline Working Directory

The pipeline requires a working directory where the FASTQ files and reference transcriptome can be accessed. 

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA

ln -s /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ/
ln -s /data/mckeeka/bulkRNA_sarcoma/reference/
```

## Generate Pipeline Configuration

The pipeline requires a working directory where the FASTQ files and reference transcriptome can be accessed. You can symlink these files instead 
of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA

ln -s /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ/
ln -s /data/mckeeka/bulkRNA_sarcoma/reference/
```






