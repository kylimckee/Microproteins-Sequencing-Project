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

### Download Human Reference Transcriptome and GTF

The human reference transcriptome (Ensembl release 115, GRCh38) was downloaded into a reference directory within bulkRNA_sarcoma.

```bash
cd /data/mckeeka/bulkRNA_sarcoma
mkdir reference
cd /data/mckeeka/bulkRNA_sarcoma/reference
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip *.gz
```

### Visualize Working Directory 

Visualize the working directory to ensure proper setup.

```bash
cd ..
ls -R
```

## Create Raw QC Pipeline Working Directory

The QC pipeline requires a working directory where the FASTQ files can be accessed. You can symlink these files instead of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir rawQC
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/rawQC
mkdir fastqc
cd /data/mckeeka/bulkRNA_sarcoma/
```

## Generate Raw QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after sequencing.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
conda create -n rawQC -c bioconda snakemake fastqc multiqc -y
conda activate rawQC
```

### Create Snakemake Raw QC Configuration File

```bash
nano rawQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("MCI_fastq_117_STS_FASTQ/{sample}.R1.fastq.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/rawQC/fastqc/{sample}.{read}_fastqc.zip", sample=SAMPLES, read=READS),
        "run_bulkRNA/rawQC/multiqc_report.html"

rule fastqc:
  input:
    "MCI_fastq_117_STS_FASTQ/{sample}.{read}.fastq.gz"
  output:
    html="run_bulkRNA/rawQC/fastqc/{sample}.{read}_fastqc.html",
    zip="run_bulkRNA/rawQC/fastqc/{sample}.{read}_fastqc.zip"
  threads: 4
  shell:
    """
    fastqc -t {threads} -o run_bulkRNA/rawQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("run_bulkRNA/rawQC/fastqc/{sample}.{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS)
  output:
    "run_bulkRNA/rawQC/multiqc_report.html"
  shell:
    """
    multiqc run_bulkRNA/rawQC/fastqc -o run_bulkRNA/rawQC
    """
```

### Run Raw QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
sbatch --cpus-per-task=4 --time=03-00:00:00 --wrap "snakemake -s rawQC_pipeline.smk -j 4"
```

## Create CutAdapt Pipeline Working Directory

The CutAdapt pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir trimmed_FASTQ
mkdir logs
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/logs
mkdir logs_CutAdapt
cd /data/mckeeka/bulkRNA_sarcoma/
```

## Generate CutAdapt Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.

### Install CutAdapt Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
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

SAMPLES = glob_wildcards("MCI_fastq_117_STS_FASTQ/{sample}.R1.fastq.gz").sample

rule all:
    input:
        expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz", sample=SAMPLES),
        expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz", sample=SAMPLES)

rule cutadapt_pe:
  input:
    r1 = "MCI_fastq_117_STS_FASTQ/{sample}.R1.fastq.gz",
    r2 = "MCI_fastq_117_STS_FASTQ/{sample}.R2.fastq.gz"
  output:
    r1 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
    r2 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
  log:
    "run_bulkRNA/logs/logs_CutAdapt/{sample}.CutAdapt.log"
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
cd /data/mckeeka/bulkRNA_sarcoma/
sbatch --cpus-per-task=4 --mem=16G --time=04-00:00:00 \--wrap "snakemake -s CutAdapt_pipeline.smk -j 4"
```

## Create Trimmed QC Pipeline Working Directory

The Trimmed QC pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir trimmedQC
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/trimmedQC
mkdir fastqc
cd /data/mckeeka/bulkRNA_sarcoma/
```

## Generate Trimmed QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after sequencing.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
conda create -n trimmedQC -c bioconda snakemake fastqc multiqc -y
conda activate trimmedQC
```

### Create Snakemake Trimmed QC Configuration File

```bash
nano trimmedQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip", sample=SAMPLES, read=READS),
        "run_bulkRNA/trimmedQC/multiqc_report.html"

rule fastqc:
  input:
    "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.{read}.trimmed.gz"
  output:
    html="run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.html",
    zip="run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip"
  threads: 4
  shell:
    """
    fastqc -t {threads} -o run_bulkRNA/trimmedQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip",
            sample=SAMPLES,
            read=READS)
  output:
    "run_bulkRNA/trimmedQC/multiqc_report.html"
  shell:
    """
    multiqc run_bulkRNA/trimmedQC/fastqc -o run_bulkRNA/trimmedQC
    """
```

### Run Trimmed QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
sbatch --cpus-per-task=4 --mem=16G --time=03-00:00:00 --wrap "snakemake -s trimmedQC_pipeline.smk -j 4"
```


## Generate Clean FASTQ Pipeline Configuration

This pipeline was generated to clean the FASTQs after sequencing to eliminate contamination.

### Install Clean FASTQ Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
conda create -n cleanFASTQ -c bioconda -c conda-forge snakemake kraken2 bowtie2 KrakenTools BEDTools samtools -y
conda activate cleanFASTQ
```

### Create Standard Databases

The Clean FASTQ pipeline requires a working directory where the standard reference databases can be accessed. You can symlink these files instead of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir clean_FASTQ
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/clean_FASTQ
mkdir kraken2_output
mkdir bowtie2_output
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/logs
mkdir logs_cleanFASTQ
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/logs/logs_cleanFASTQ
mkdir logs_kraken2
mkdir logs_bowtie2

cd /data/mckeeka/bulkRNA_sarcoma/reference

#Copy Kraken2 Standard Reference Database from Biowulf
cp -r /fdb/kraken/20220803_standard_kraken2 kraken2_database

#Create Contaminant Reference for Bowtie2
grep -E 'gene_biotype "(artifact|Mt_rRNA|Mt_tRNA|ribozyme|rRNA|rRNA_pseudogene|scaRNA|snoRNA|snRNA|vault_RNA)"' Homo_sapiens.GRCh38.115.gtf > contaminants.gtf

awk '$3=="exon"' contaminants.gtf | \
awk '{print$1"\t"$4-1"\t"$5}' > contaminants.bed

bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed contaminants.bed -fo contaminants.fa

mkdir -p contaminants_index
bowtie2-build contaminants.fa contaminants_index/contaminants
```

### Create Snakemake Clean FASTQ Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma
nano cleanFASTQ_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/clean_FASTQ/{sample}.fastq.{read}.clean.gz", sample=SAMPLES, read=READS)

rule kraken2:
  input:
    r1 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
    r2 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
  output:
    report = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_report.txt",
    output = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_output.txt"
  threads: 4
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_kraken2/{sample}.kraken2.log"
  shell:
    """
    kraken2 \
        --db reference/kraken2_database \
        --paired {input.r1} {input.r2} \
        --report {output.report} \
        --output {output.output} \
        --threads {threads} \
        &> {log}
    """

rule extract_human_unclassified:
    input:
        kraken2="run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_output.txt",
        report="run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_report.txt",
        r1="run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
        r2="run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
    output:
        r1="run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R1.kraken.gz",
        r2="run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R2.kraken.gz"
    log:
        "run_bulkRNA/logs/logs_cleanFASTQ/logs_kraken2/{sample}.kraken2_filter.log"
    shell:
        r"""
        set -euo pipefail

        tmp_h1=$(mktemp)
        tmp_h2=$(mktemp)
        tmp_u1=$(mktemp)
        tmp_u2=$(mktemp)

        # Extract human reads
        extract_kraken_reads.py \
            -k {input.kraken2} \
            -r {input.report} \
            -s {input.r1} \
            -s2 {input.r2} \
            -t 9606 \
            --include-children \
            --fastq-output \
            -o "$tmp_h1" \
            -o2 "$tmp_h2"

        # Extract unclassified reads
        extract_kraken_reads.py \
            -k {input.kraken2} \
            -r {input.report} \
            -s {input.r1} \
            -s2 {input.r2} \
            -t 0 \
            --fastq-output \
            -o "$tmp_u1" \
            -o2 "$tmp_u2"

        # Ensure files exist (handles edge cases)
        [ -s "$tmp_h1" ] || touch "$tmp_h1"
        [ -s "$tmp_h2" ] || touch "$tmp_h2"
        [ -s "$tmp_u1" ] || touch "$tmp_u1"
        [ -s "$tmp_u2" ] || touch "$tmp_u2"

        # Combine
        cat "$tmp_h1" "$tmp_u1" | gzip > {output.r1}
        cat "$tmp_h2" "$tmp_u2" | gzip > {output.r2}

        rm -f "$tmp_h1" "$tmp_h2" "$tmp_u1" "$tmp_u2"
        """ + " &> {log}"

rule bowtie2_contaminant_mapping:
  input:
    r1 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R1.kraken.gz",
    r2 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R2.kraken.gz"
  output:
    bam = "run_bulkRNA/clean_FASTQ/bowtie2_output/{sample}_contamination.bam"
  threads: 4
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_bowtie2/{sample}.bowtie2.log"
  shell:
    """
    bowtie2 \
        -x reference/contaminants_index/contaminants \
        -1 {input.r1} \
        -2 {input.r2} \
        --sensitive \
        --threads {threads} \
        2> {log} \
        | samtools view -b -o {output.bam} -
    """

rule filter_unmapped:
  input:
    bam = "run_bulkRNA/clean_FASTQ/bowtie2_output/{sample}_contamination.bam"
  output:
    r1 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R1.clean.gz",
    r2 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R2.clean.gz"
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_bowtie2/{sample}.bowtie2_filter.log"
  shell:
    r"""
    set -euo pipefail

    tmp_bam=$(mktemp --suffix=.bam)
    tmp_namesort=$(mktemp --suffix=.bam)
    tmp_r1=$(mktemp --suffix=.fq)
    tmp_r2=$(mktemp --suffix=.fq)

    samtools view -b -f 12 -F 256 {input.bam} > "$tmp_bam"
    samtools sort -n -o "$tmp_namesort" "$tmp_bam"

    bedtools bamtofastq \
        -i "$tmp_namesort" \
        -fq "$tmp_r1" \
        -fq2 "$tmp_r2"

    gzip -c "$tmp_r1" > {output.r1}
    gzip -c "$tmp_r2" > {output.r2}

    rm -f "$tmp_bam" "$tmp_namesort" "$tmp_r1" "$tmp_r2"
    """ + "&> {log}"

```

### Run Clean FASTQ Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma
sbatch --cpus-per-task=4 --mem=64G --time=10-00:00:00 \--wrap "snakemake -s cleanFASTQ_pipeline.smk -j 4"
```


## Create Read Counts QC Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/clean_FASTQ/
mkdir ReadCounts_output
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/
```


## Generate Read Counts Pipeline Configuration

This pipeline was generated to count the reads that were filtered out of each step of the Clean FASTQ pipeline.

### Install Read Counts Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
conda create -n ReadCounts -c conda-forge -c bioconda snakemake python=3.10 pandas -y
conda activate ReadCounts
```

### Create Snakemake Read Counts Configuration File

```bash
nano ReadCounts_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("trimmed_FASTQ/{sample}.fastq.{read}.trimmed.gz").sample

rule all:
    input:
        "clean_FASTQ/ReadCounts_output/read_summary.csv",
        expand("trimmed_FASTQ/{sample}.fastq.{read}.count.txt", sample=SAMPLES, read=[1,2]),
        expand("clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.count.txt", sample=SAMPLES, read=[1,2]),
        expand("clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.count.txt", sample=SAMPLES, read=[1,2])

rule count_trimmed_reads:
  input:
    "trimmed_FASTQ/{sample}.fastq.{read}.trimmed.gz"
  output:
    "trimmed_FASTQ/{sample}.fastq.{read}.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule count_kraken_reads:
  input:
    "clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.gz"
  output:
    "clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule count_bowtie_reads:
  input:
    "clean_FASTQ/{sample}.fastq.{read}.clean.gz"
  output:
    "clean_FASTQ/{sample}.fastq.{read}.clean.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule summarize_read_counts:
  input:
      TRIMMED = expand("trimmed_FASTQ/{sample}.fastq.{read}.count.txt", sample=SAMPLES, read=[1,2]),
      KRAKEN = expand("clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.count.txt", sample=SAMPLES, read=[1,2]),
      BOWTIE = expand("clean_FASTQ/{sample}.fastq.{read}.clean.count.txt", sample=SAMPLES, read=[1,2])
  output:
      "clean_FASTQ/ReadCounts_output/read_summary.csv"
  run:
      import pandas a pd
      data = []
      for sample in SAMPLES:
        row = {"sample": sample}
        for step, path_template in [
            ("TRIMMED", "trimmed_FASTQ/{sample}.fastq.{read}.count.txt"),
            ("KRAKEN", "clean_FASTQ/kraken2_output/{sample}.fastq.{read}.kraken.count.txt"),
            ("BOWTIE", "clean_FASTQ/{sample}.fastq.{read}.clean.count.txt")
        ]:
            for read in [1,2]:
              col_name = f"{step}_{read}"
              file_path = path_template.format(sample=sample, read=read)
              row[col_name] = int(float(open(file_path).read().strip()))
        data.append(row)
      df = pd.DataFrame(data)
      df.to_csv(output[0], index=False)

```

### Run Read Counts Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
sbatch --cpus-per-task=4 --mem=16G --time=04:00:00 \--wrap "snakemake -s ReadCounts_pipeline.smk -j 4"
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

## Create STAR Mapping Pipeline Working Directory

The STAR Mapping pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA
mkdir STAR
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/logs
mkdir logs_STAR
cd /data/mckeeka/bulkRNA_sarcoma/
```

## Generate STAR Mapping Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.

### Install STAR Mapping Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n STARmap -c bioconda snakemake star SAMtools -y
conda activate STARmap
```

### Create Snakemake STAR Mapping Configuration File

```bash
nano STARmap_pipeline.smk

# Add the following code to the configuration file:

GENOME_FASTA = "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
STAR_INDEX_DIR = "reference/STAR_index"

SAMPLES = glob_wildcards("MCI_fastq_117_STS_FASTQ/{sample}.R1.fastq.gz").sample

rule all:
  input:
    "reference/STAR_index",
    expand("run_bulkRNA/STAR/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

rule star_index:
  input:
    fasta=GENOME_FASTA,
    gtf=GTF
  output:
    STAR_INDEX_DIR
  params:
    outdir=STAR_INDEX_DIR
  threads: 4
  shell:
    """
    mkdir -p {params.outdir}

    STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99
    """

rule star_two_pass:
  input:
    r1 = "MCI_fastq_117_STS_FASTQ/{sample}.R1.fastq.gz",
    r2 = "MCI_fastq_117_STS_FASTQ/{sample}.R2.fastq.gz",
    index = STAR_INDEX_DIR
  output:
    bam = "run_bulkRNA/STAR/{sample}.Aligned.sortedByCoord.out.bam"
  log:
    "run_bulkRNA/logs/logs_STAR/{sample}.log"
  threads: 4
  shell:
      """
      STAR \
          --runThreadN {threads} \
          --genomeDir {input.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --twopassMode Basic \
          --chimSegmentMin 12 \
          --outFilterMultimapNmax 20 \
          --winAnchorMultimapNmax 50 \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix run_bulkRNA/STAR/{wildcards.sample}. \
          &> {log}
      """

rule samtools_index:
  input:
    bam = "run_bulkRNA/STAR/{sample}.Aligned.sortedByCoord.out.bam"
  output:
    bai = "run_bulkRNA/STAR/{sample}.Aligned.sortedByCoord.out.bam.bai"
  threads: 4
  shell:
    """
    samtools index {input.bam}
    """

```

### Run STAR Mapping Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma/
sbatch --cpus-per-task=4 --mem=64G --time=10-00:00:00 --wrap "snakemake -s STARmap_pipeline.smk --cores 4"
```

## Generate STAR QC Pipeline Configuration

This pipeline was generated to analyze the quality of STAR mapping on our unfiltered samples.

### Install STAR QC Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma
conda create -n STAR_QC -c conda-forge -c bioconda snakemake python=3.10 -y
conda activate STAR_QC
```

### Create Snakemake STAR QC Configuration File

```bash
nano STAR_QC.py


# Add the following code to the configuration file:

import csv
import glob
import os

STAR_LOG_DIR = "run_bulkRNA/STAR"
OUT_CSV = "run_bulkRNA/STAR/star_qc_summary.csv"

def parse_log_final_out(path):
    # Example filename: PBCKHM_0E7M4E.Log.final.out -> PBCKHM_0E7M4E
    sample = os.path.basename(path).replace(".Log.final.out", "")
    metrics = {"sample": sample}

    with open(path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = line.split("|", 1)
            key = key.strip()
            val = val.strip()
            metrics[key] = val

    return metrics

logs = glob.glob(os.path.join(STAR_LOG_DIR, "*.Log.final.out"))
rows = [parse_log_final_out(p) for p in logs]

if not rows:
    raise SystemExit(f"No STAR log files found in {STAR_LOG_DIR}")

all_keys = sorted({k for row in rows for k in row.keys()})
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_keys)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {OUT_CSV}")

```

### Run STAR QC Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma
python STAR_QC.py
```
## Create FeatureCounts Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/
mkdir FeatureCounts
mkdir TPM
```

## Generate Feature Counts Pipeline Configuration

This pipeline was generated to count the features after STAR two-pass mapping on the unfiltered samples.

### Install Feature Counts Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma
conda create -n FeatureCounts -c conda-forge -c bioconda snakemake subread python=3.10 -y
conda activate FeatureCounts
```

### Create Snakemake Feature Counts Configuration File

```bash
nano FeatureCounts_pipeline.smk

# Add the following code to the configuration file:

from pathlib import Path
from glob import glob

GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
BAM_DIR = "run_bulkRNA/STAR"
COUNT_DIR = "run_bulkRNA/FeatureCounts"
TPM_DIR = "run_bulkRNA/TPM"

SAMPLES = sorted(
    Path(b).name.replace(".Aligned.sortedByCoord.out.bam", "")
    for b in glob(f"{BAM_DIR}/*.Aligned.sortedByCoord.out.bam")
)

rule all:
    input:
        f"{COUNT_DIR}/gene_counts_unfiltered.txt",
        f"{COUNT_DIR}/gene_counts_clean_unfiltered.txt",
        f"{TPM_DIR}/gene_tpm_unfiltered.tsv"

rule featurecounts:
    input:
        bams=lambda wildcards: expand(
            f"{BAM_DIR}" + "/{sample}.Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        ),
        gtf=GTF
    output:
        counts=f"{COUNT_DIR}/gene_counts_unfiltered.txt"
    threads: 6
    shell:
        """
        mkdir -p {COUNT_DIR}
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -p \
            -B \
            -C \
            {input.bams}
        """

rule clean_featurecounts:
    input:
        counts=f"{COUNT_DIR}/gene_counts_unfiltered.txt"
    output:
        clean=f"{COUNT_DIR}/gene_counts_clean_unfiltered.txt"
    run:
        import pandas as pd
        import re

        df = pd.read_csv(
            input.counts,
            sep="\t",
            comment="#"
        )

        annotation_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        sample_cols = [c for c in df.columns if c not in annotation_cols]

        clean_names = [re.sub(r".*/", "", c) for c in sample_cols]
        clean_names = [re.sub(r"\.Aligned\.sortedByCoord\.out\.bam$", "", c) for c in clean_names]

        rename_map = dict(zip(sample_cols, clean_names))
        df = df.rename(columns=rename_map)

        df.to_csv(output.clean, sep="\t", index=False)

rule tpm:
    input:
        counts=f"{COUNT_DIR}/gene_counts_clean_unfiltered.txt"
    output:
        tpm=f"{TPM_DIR}/gene_tpm_unfiltered.tsv"
    run:
        import pandas as pd
        import os

        os.makedirs(TPM_DIR, exist_ok=True)

        df = pd.read_csv(
            input.counts,
            sep="\t",
            comment="#"
        )

        annotation_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        sample_cols = [c for c in df.columns if c not in annotation_cols]

        lengths_kb = df["Length"] / 1000.0
        counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

        rpk = counts.div(lengths_kb, axis=0)
        scaling_factors = rpk.sum(axis=0)
        tpm = rpk.div(scaling_factors, axis=1) * 1e6

        tpm.insert(0, "Geneid", df["Geneid"])
        tpm.to_csv(output.tpm, sep="\t", index=False)
```

### Run Feature Counts Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_sarcoma
sbatch --cpus-per-task=6 --time=02:00:00 --wrap "snakemake -s FeatureCounts_pipeline.smk --cores 6"
```

## Create Microproteins Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/
mkdir Microproteins
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/Microproteins
mkdir high_confidence
mkdir discovery
```

## Generate Microproteins Pipeline Configuration

### Install Microproteins Tools

```bash
cd /data/mckeeka/bulkRNA_sarcoma
conda create -n Microproteins -c conda-forge -c bioconda snakemake gffread transdecoder stringtie r-base python=3.10 -y
conda activate Microproteins
```

### Create High Confidence Microproteins Configuration File

```bash
nano Microproteins_highconfidence.smk

# Add the following code to the configuration file:

from pathlib import Path
from glob import glob

MIN_MICROPROTEIN_AA = 100
GENOME_FASTA = "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
COUNTS = "run_bulkRNA/FeatureCounts/gene_counts_filtered.txt"
OUT_DIR = "run_bulkRNA/Microproteins/high_confidence"
BAM_DIR = "run_bulkRNA/STAR"

SAMPLES = sorted(
    Path(b).name.replace(".Aligned.sortedByCoord.out.bam", "")
    for b in glob(f"{BAM_DIR}/*.Aligned.sortedByCoord.out.bam")
)

rule all:
    input:
        expand("run_bulkRNA/Microproteins/stringtie/{sample}.gtf", sample=SAMPLES),
        "run_bulkRNA/Microproteins/high_confidence/orf/merged.gtf",
        f"{OUT_DIR}/orf/transcripts.fa",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.gff3",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.cds",
        f"{OUT_DIR}/orf/orf_metadata.tsv",
        f"{OUT_DIR}/orf/microproteins.tsv"

rule stringtie_assemble:
    input:
        bam=f"{BAM_DIR}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        gtf="run_bulkRNA/Microproteins/stringtie/{sample}.gtf"
    shell:
        """
        mkdir -p run_bulkRNA/Microproteins/stringtie
        stringtie {input.bam} -G {GTF} -o {output.gtf}
        """

rule merge_transcripts:
    input:
        gtfs=expand("run_bulkRNA/Microproteins/stringtie/{sample}.gtf", sample=SAMPLES)
    output:
        merged="run_bulkRNA/Microproteins/high_confidence/orf/merged.gtf"
    shell:
        """
        mkdir -p run_bulkRNA/Microproteins/high_confidence/orf
        stringtie --merge -G {GTF} -o {output.merged} {input.gtfs}
        """

rule extract_transcripts:
    input:
        genome=GENOME_FASTA,
        gtf="run_bulkRNA/Microproteins/high_confidence/orf/merged.gtf"
    output:
        fa=f"{OUT_DIR}/orf/transcripts.fa"
    shell:
        """
        mkdir -p {OUT_DIR}/orf
        gffread {input.gtf} -g {input.genome} -w {output.fa}
        """

rule transdecoder_longorfs:
    input:
        fa=f"{OUT_DIR}/orf/transcripts.fa"
    output:
        flag=f"{OUT_DIR}/orf/.longorfs.done"
    shell:
        """
        cd {OUT_DIR}/orf
        TransDecoder.LongOrfs -t transcripts.fa
        touch .longorfs.done
        """

rule transdecoder_predict:
    input:
        fa=f"{OUT_DIR}/orf/transcripts.fa",
        flag=f"{OUT_DIR}/orf/.longorfs.done"
    output:
        pep=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep",
        gff3=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.gff3",
        cds=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.cds"
    shell:
        """
        cd {OUT_DIR}/orf
        TransDecoder.Predict -t transcripts.fa
        """

rule build_orf_metadata:
    input:
        pep=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep"
    output:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    run:
        import re

        def parse_fasta(path):
            header = None
            seq_parts = []
            with open(path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq_parts)
                        header = line[1:]
                        seq_parts = []
                    else:
                        seq_parts.append(line)
                if header is not None:
                    yield header, "".join(seq_parts)

        with open(output.tsv, "w") as out:
            out.write("orf_id\ttranscript_id\tlength_aa\tsequence\tcategory\n")
            for header, seq in parse_fasta(input.pep):
                orf_id = header.split()[0]

                transcript_id = orf_id
                if ".p" in orf_id:
                    transcript_id = orf_id.rsplit(".p", 1)[0]

                m = re.search(r"len=(\d+)", header)
                length_aa = int(m.group(1)) if m else len(seq)

                category = "microprotein" if length_aa <= MIN_MICROPROTEIN_AA else "protein"

                out.write(f"{orf_id}\t{transcript_id}\t{length_aa}\t{seq}\t{category}\n")


rule split_microproteins:
    input:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    output:
        tsv=f"{OUT_DIR}/orf/microproteins.tsv"
    run:
        import csv

        with open(input.tsv) as inf, open(output.tsv, "w", newline="") as outf:
            reader = csv.DictReader(inf, delimiter="\t")
            writer = csv.DictWriter(outf, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if row["category"] == "microprotein":
                    writer.writerow(row)


rule count_qc_plots:
    input:
        counts=COUNTS
    output:
        lib="run_bulkRNA/Microproteins/qc/library_sizes.pdf",
        dist="run_bulkRNA/Microproteins/qc/count_distribution.pdf"
    shell:
        r"""
        mkdir -p run_bulkRNA/Microproteins/qc
        Rscript -e '
          counts <- read.delim(
            "{input.counts}",
            header = TRUE,
            sep = "\t",
            comment.char = "#",
            check.names = FALSE,
            fill = TRUE,
            quote = "",
            stringsAsFactors = FALSE
          )

          annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
          sample_counts <- counts[, !(names(counts) %in% annotation_cols), drop = FALSE]

          if ("Geneid" %in% names(counts)) {{
            rownames(sample_counts) <- counts$Geneid
          }}

          sample_counts <- as.data.frame(lapply(sample_counts, as.numeric))

          pdf("{output.lib}")
          barplot(colSums(sample_counts, na.rm = TRUE),
                  las = 2,
                  main = "Library sizes",
                  ylab = "Total counts")
          dev.off()

          pdf("{output.dist}")
          boxplot(log2(sample_counts + 1),
                  las = 2,
                  main = "Log2 count distribution",
                  ylab = "log2(counts + 1)")
          dev.off()
        '
        """
```

### Run Microproteins High Confidence Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
sbatch --time=00-08:00:00  --cpus-per-task=8 --mem=32G --wrap="snakemake -s Microproteins_highconfidence.smk --cores 8"
```

### Create Discovery Microproteins Configuration File

```bash
nano Microproteins_discovery.smk

# Add the following code to the configuration file:

from pathlib import Path

MIN_MICROPROTEIN_AA = 100
pep_fasta = "run_bulkRNA/Microproteins/high_confidence/orf/transcripts.fa.transdecoder_dir/longest_orfs.pep"
OUT_DIR = "run_bulkRNA/Microproteins/discovery"

rule all:
    input:
        f"{OUT_DIR}/orf/orf_metadata.tsv",
        f"{OUT_DIR}/orf/microproteins.tsv"

rule build_orf_metadata:
    input:
        pep=pep_fasta
    output:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    run:
        import re

        Path(f"{OUT_DIR}/orf").mkdir(parents=True, exist_ok=True)

        def parse_fasta(path):
            header = None
            seq_parts = []
            with open(path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq_parts)
                        header = line[1:]
                        seq_parts = []
                    else:
                        seq_parts.append(line)
                if header is not None:
                    yield header, "".join(seq_parts)

        with open(output.tsv, "w") as out:
            out.write("orf_id\ttranscript_id\tlength_aa\tsequence\tcategory\n")
            for header, seq in parse_fasta(input.pep):
                orf_id = header.split()[0]

                transcript_id = orf_id
                if ".p" in orf_id:
                    transcript_id = orf_id.rsplit(".p", 1)[0]

                m = re.search(r"len=(\d+)", header)
                length_aa = int(m.group(1)) if m else len(seq)
                category = "microprotein" if length_aa <= MIN_MICROPROTEIN_AA else "protein"

                out.write(f"{orf_id}\t{transcript_id}\t{length_aa}\t{seq}\t{category}\n")

rule split_microproteins:
    input:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    output:
        tsv=f"{OUT_DIR}/orf/microproteins.tsv"
    run:
        import csv

        with open(input.tsv) as inf, open(output.tsv, "w", newline="") as outf:
            reader = csv.DictReader(inf, delimiter="\t")
            writer = csv.DictWriter(outf, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if row["category"] == "microprotein":
                    writer.writerow(row)
```

### Run Discovery Microproteins Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
sbatch --time=00-00:20:00  --cpus-per-task=8 --mem=32G --wrap="snakemake -s Microproteins_discovery.smk --cores 8"
```

## Create Novel Microproteins Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/
mkdir Novel_Microproteins
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA/Novel_Microproteins/
mkdir high_confidence
mkdir discovery
```

## Generate Novel Microproteins Pipeline Configuration

### Install Novel Microproteins Tools

```bash
conda create -n smorf \
  -c conda-forge -c bioconda \
  r-base=4.4 \
  r-data.table \
  r-stringr \
  bioconductor-biostrings \
  bioconductor-rtracklayer \
  diamond \
  -y

conda activate smorf
```

### Create Novel Microproteins High Confidence Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma

nano Microproteins_Novel_HighConfidence.r

# Add the following code to the configuration file:

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(stringr)
  library(rtracklayer)
})

# =========================
# User inputs
# =========================

pep_fasta    <- "run_bulkRNA/Microproteins/high_confidence/orf/transcripts.fa.transdecoder.pep"
gtf_file     <- "reference/Homo_sapiens.GRCh38.115.gtf"

out_dir <- "run_bulkRNA/Novel_Microproteins/high_confidence"
out_prefix   <- file.path(out_dir, "smorf_results")
max_aa_len   <- 100

diamond_evalue    <- "1e-5"
diamond_max_hits   <- 5

# DIAMOND databases
uniprot_full_db <- "reference/uniprot_full.dmnd"
sprot_db        <- "reference/uniprot_sprot.dmnd"
openprot_db     <- "reference/openprot_human.dmnd"
smprot_db       <- "reference/smprot2_human_ribo.dmnd"

# =========================
# Helper functions
# =========================

parse_transdecoder_header <- function(hdr) {
  # Example:
  # ENST00000832824.p1 type:complete gc:universal ENST00000832824:764-345(-)

  first_token <- sub(" .*", "", hdr)

  tx_orf <- first_token
  tx_id <- sub("\\.p[0-9]+$", "", tx_orf)
  orf_id <- tx_orf

  coord_str <- NA_character_
  if (grepl(":[0-9]+-[0-9]+\\([+-]\\)$", hdr)) {
    coord_str <- sub(".* ([^ ]+:[0-9]+-[0-9]+\\([+-]\\))$", "\\1", hdr)
  }

  transcript_id_from_coord <- NA_character_
  tx_start <- NA_integer_
  tx_end <- NA_integer_
  strand <- NA_character_

  if (!is.na(coord_str)) {
    m <- str_match(coord_str, "^(.+):([0-9]+)-([0-9]+)\\(([+-])\\)$")
    transcript_id_from_coord <- m[, 2]
    tx_start <- suppressWarnings(as.integer(m[, 3]))
    tx_end <- suppressWarnings(as.integer(m[, 4]))
    strand <- m[, 5]
  }

  data.table(
    header = hdr,
    orf_id = orf_id,
    tx_id = tx_id,
    tx_id_from_coord = transcript_id_from_coord,
    tx_start = tx_start,
    tx_end = tx_end,
    strand = strand
  )
}

clean_interval <- function(a, b) {
  c(min(a, b), max(a, b))
}

map_tx_interval_to_genome <- function(tx_start, tx_end, exons_df) {
  if (nrow(exons_df) == 0 || any(is.na(c(tx_start, tx_end)))) {
    return(list(
      chrom = NA_character_,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = NA_character_
    ))
  }

  exons_df <- as.data.table(exons_df)

  strand <- unique(exons_df$strand)
  chrom <- unique(as.character(exons_df$seqnames))
  strand <- strand[1]
  chrom <- chrom[1]

  if (strand == "+") {
    setorder(exons_df, start, end)
  } else {
    setorder(exons_df, -start, -end)
  }

  exons_df[, exon_len := end - start + 1L]
  exons_df[, tx_exon_start := cumsum(c(0L, head(exon_len, -1L))) + 1L]
  exons_df[, tx_exon_end := cumsum(exon_len)]

  map_one_pos <- function(tx_pos) {
    hit <- exons_df[tx_exon_start <= tx_pos & tx_exon_end >= tx_pos][1]
    if (nrow(hit) == 0) return(NA_integer_)

    offset <- tx_pos - hit$tx_exon_start
    if (strand == "+") {
      hit$start + offset
    } else {
      hit$end - offset
    }
  }

  g1 <- map_one_pos(tx_start)
  g2 <- map_one_pos(tx_end)

  if (is.na(g1) || is.na(g2)) {
    return(list(
      chrom = chrom,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = strand
    ))
  }

  list(
    chrom = chrom,
    genomic_start = min(g1, g2),
    genomic_end = max(g1, g2),
    strand = strand
  )
}

run_diamond <- function(query_fa, db, out_tsv, evalue = "1e-5", max_hits = 5) {
  diamond_cmd <- c(
    "blastp",
    "-q", query_fa,
    "-d", db,
    "-o", out_tsv,
    "-e", evalue,
    "-k", as.character(max_hits),
    "--quiet",
    "--outfmt", "6",
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )

  message("Running DIAMOND against: ", db)
  status <- system2("diamond", diamond_cmd)
  if (!identical(status, 0L)) {
    warning("DIAMOND returned non-zero exit status for: ", db)
  }
}

read_best_hits <- function(tsv_file, db_name) {
  if (!file.exists(tsv_file) || file.info(tsv_file)$size == 0) return(NULL)

  hits <- fread(tsv_file, header = FALSE)
  setnames(hits, c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  ))

  hits[, `:=`(
    pident = as.numeric(pident),
    length = as.numeric(length),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    qlen = as.numeric(qlen),
    slen = as.numeric(slen)
  )]

  hits[, `:=`(
    qcov = length / qlen,
    scov = length / slen,
    db = db_name
  )]

  # Best hit = lowest evalue, then highest bitscore, then highest identity, then highest qcov
  setorder(hits, qseqid, evalue, -bitscore, -pident, -qcov)
  hits[, .SD[1], by = qseqid]
}

# =========================
# Read and filter peptides
# =========================

pep <- readAAStringSet(pep_fasta)
hdrs <- names(pep)
aa_len <- width(pep)

meta <- rbindlist(lapply(hdrs, parse_transdecoder_header), fill = TRUE)
meta[, aa_len := aa_len]

# Keep only ORFs below the threshold
keep_idx <- aa_len <= max_aa_len
pep_smorf <- pep[keep_idx]
meta_smorf <- meta[keep_idx]

# Write filtered peptide FASTA
smorf_pep_fa <- paste0(out_prefix, ".smorfs_lt", max_aa_len, "aa.pep.fa")
names(pep_smorf) <- meta_smorf$orf_id
writeXStringSet(pep_smorf, smorf_pep_fa)

# =========================
# DIAMOND searches
# =========================

openprot_tsv <- paste0(out_prefix, ".openprot.tsv")
smprot_tsv   <- paste0(out_prefix, ".smprot.tsv")
sprot_tsv    <- paste0(out_prefix, ".uniprot_sprot.tsv")
uniprot_tsv  <- paste0(out_prefix, ".uniprot_full.tsv")

run_diamond(smorf_pep_fa, openprot_db, openprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, smprot_db,   smprot_tsv,   evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, sprot_db,    sprot_tsv,    evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, uniprot_full_db, uniprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)

openprot_best <- read_best_hits(openprot_tsv, "openprot")
smprot_best   <- read_best_hits(smprot_tsv,   "smprot")
sprot_best    <- read_best_hits(sprot_tsv,    "sprot")
uniprot_best  <- read_best_hits(uniprot_tsv,  "uniprot_full")

# =========================
# GTF processing
# =========================

tx_anno <- NULL
tx_exons <- NULL

if (!is.null(gtf_file) && file.exists(gtf_file)) {
  message("Reading GTF...")
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  rm(gtf)
  gc()

  if (!("transcript_id" %in% names(gtf_df))) {
    stop("GTF does not contain a transcript_id column in metadata.")
  }

  anno_cols <- intersect(
    c("transcript_id", "gene_id", "gene_name", "gene_type", "transcript_type", "biotype"),
    names(gtf_df)
  )

  if (length(anno_cols) > 0) {
    tx_anno <- unique(gtf_df[, anno_cols, drop = FALSE])
  }

  exons_only <- gtf_df[gtf_df$type == "exon", , drop = FALSE]

  if (nrow(exons_only) > 0) {
    exons_only <- exons_only[, intersect(
      c("seqnames", "start", "end", "strand", "transcript_id"),
      names(exons_only)
    ), drop = FALSE]

    tx_exons <- as.data.table(exons_only)
    setkey(tx_exons, transcript_id)
  }
}

# =========================
# Build annotation table
# =========================

smorf_dt <- copy(meta_smorf)

# Use transcript ID from coordinate field if available
smorf_dt[, transcript_id := fifelse(!is.na(tx_id_from_coord), tx_id_from_coord, tx_id)]

# Initialize columns
smorf_dt[, `:=`(
  openprot_hit = NA_character_,
  openprot_pident = NA_real_,
  openprot_qcov = NA_real_,
  openprot_evalue = NA_real_,
  openprot_bitscore = NA_real_,

  smprot_hit = NA_character_,
  smprot_pident = NA_real_,
  smprot_qcov = NA_real_,
  smprot_evalue = NA_real_,
  smprot_bitscore = NA_real_,

  sprot_hit = NA_character_,
  sprot_pident = NA_real_,
  sprot_qcov = NA_real_,
  sprot_evalue = NA_real_,
  sprot_bitscore = NA_real_,

  uniprot_hit = NA_character_,
  uniprot_pident = NA_real_,
  uniprot_qcov = NA_real_,
  uniprot_evalue = NA_real_,
  uniprot_bitscore = NA_real_,

  genomic_chrom = NA_character_,
  genomic_start = NA_integer_,
  genomic_end = NA_integer_,
  genomic_strand = NA_character_
)]

# Attach DIAMOND hits by qseqid -> orf_id
if (!is.null(openprot_best) && nrow(openprot_best) > 0) {
  smorf_dt[openprot_best, on = c(orf_id = "qseqid"), `:=`(
    openprot_hit = i.sseqid,
    openprot_pident = i.pident,
    openprot_qcov = i.qcov,
    openprot_evalue = i.evalue,
    openprot_bitscore = i.bitscore
  )]
}

if (!is.null(smprot_best) && nrow(smprot_best) > 0) {
  smorf_dt[smprot_best, on = c(orf_id = "qseqid"), `:=`(
    smprot_hit = i.sseqid,
    smprot_pident = i.pident,
    smprot_qcov = i.qcov,
    smprot_evalue = i.evalue,
    smprot_bitscore = i.bitscore
  )]
}

if (!is.null(sprot_best) && nrow(sprot_best) > 0) {
  smorf_dt[sprot_best, on = c(orf_id = "qseqid"), `:=`(
    sprot_hit = i.sseqid,
    sprot_pident = i.pident,
    sprot_qcov = i.qcov,
    sprot_evalue = i.evalue,
    sprot_bitscore = i.bitscore
  )]
}

if (!is.null(uniprot_best) && nrow(uniprot_best) > 0) {
  smorf_dt[uniprot_best, on = c(orf_id = "qseqid"), `:=`(
    uniprot_hit = i.sseqid,
    uniprot_pident = i.pident,
    uniprot_qcov = i.qcov,
    uniprot_evalue = i.evalue,
    uniprot_bitscore = i.bitscore
  )]
}

# =========================
# Optional genomic mapping
# =========================

if (!is.null(tx_exons)) {
  map_res <- lapply(seq_len(nrow(smorf_dt)), function(i) {
    tx <- smorf_dt$transcript_id[i]
    if (is.na(tx) || !(tx %in% tx_exons$transcript_id)) {
      return(list(
        chrom = NA_character_,
        genomic_start = NA_integer_,
        genomic_end = NA_integer_,
        strand = NA_character_
      ))
    }

    tx_df <- tx_exons[J(tx), nomatch = 0L]
    tx_df <- as.data.frame(tx_df)
    tx_df <- tx_df[, intersect(c("seqnames", "start", "end", "strand"), names(tx_df)), drop = FALSE]

    interval <- clean_interval(smorf_dt$tx_start[i], smorf_dt$tx_end[i])
    map_tx_interval_to_genome(interval[1], interval[2], tx_df)
  })

smorf_dt[, genomic_chrom := vapply(map_res, function(x) as.character(x$chrom), character(1))]
smorf_dt[, genomic_start := vapply(map_res, function(x) as.integer(x$genomic_start), integer(1))]
smorf_dt[, genomic_end := vapply(map_res, function(x) as.integer(x$genomic_end), integer(1))]
smorf_dt[, genomic_strand := vapply(map_res, function(x) as.character(x$strand), character(1))]
}

# =========================
# Classification logic
# =========================

smorf_dt[, `:=`(
  hit_in_openprot = !is.na(openprot_hit),
  hit_in_smprot = !is.na(smprot_hit),

  strong_hit_in_sprot = !is.na(sprot_hit) &
    !is.na(sprot_evalue) & sprot_evalue <= 1e-20 &
    !is.na(sprot_pident) & sprot_pident >= 90 &
    !is.na(sprot_qcov) & sprot_qcov >= 0.80,

  weak_hit_in_uniprot = !is.na(uniprot_hit) &
    !is.na(uniprot_evalue) & uniprot_evalue <= 1e-5 &
    !is.na(uniprot_pident) & uniprot_pident >= 30 &
    !is.na(uniprot_qcov) & uniprot_qcov >= 0.50
)]

smorf_dt[, class := fifelse(
  hit_in_openprot | hit_in_smprot,
  "known_microprotein",
  fifelse(
    strong_hit_in_sprot,
    "known_canonical",
    fifelse(
      weak_hit_in_uniprot,
      "homologous",
      "novel_candidate"
    )
  )
)]

# =========================
# Final cleanup / output
# =========================

# Optional: order columns a little more nicely
setcolorder(smorf_dt, c(
  "header", "orf_id", "tx_id", "tx_id_from_coord", "transcript_id",
  "tx_start", "tx_end", "strand", "aa_len",
  "genomic_chrom", "genomic_start", "genomic_end", "genomic_strand",
  "openprot_hit", "openprot_pident", "openprot_qcov", "openprot_evalue", "openprot_bitscore",
  "smprot_hit", "smprot_pident", "smprot_qcov", "smprot_evalue", "smprot_bitscore",
  "sprot_hit", "sprot_pident", "sprot_qcov", "sprot_evalue", "sprot_bitscore",
  "uniprot_hit", "uniprot_pident", "uniprot_qcov", "uniprot_evalue", "uniprot_bitscore",
  "hit_in_openprot", "hit_in_smprot", "strong_hit_in_sprot", "weak_hit_in_uniprot",
  "class"
))

fwrite(
  smorf_dt,
  paste0(out_prefix, ".smorf_classification.tsv"),
  sep = "\t",
  na = "NA"
)

message("Done. Wrote: ", paste0(out_prefix, ".smorf_classification.tsv"))
```

### Run Novel Microproteins High Confidence Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma/

sbatch --job-name=smorf --cpus-per-task=4 --mem=32G --time=4:00:00 --output=smorf_%j.out --error=smorf_%j.err --wrap="set -x; echo START $(date); source ~/.bashrc; conda activate smorf; which Rscript; echo ENV_OK; Rscript Microproteins_Novel_HighConfidence.r; echo DONE $(date)"

```

### Create Novel Microproteins Discovery Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma

nano Microproteins_Novel_Discovery.r

# Add the following code to the configuration file:

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(stringr)
  library(rtracklayer)
})

# =========================
# User inputs
# =========================

pep_fasta    <- "run_bulkRNA/Microproteins/high_confidence/orf/transcripts.fa.transdecoder_dir/longest_orfs.pep"
gtf_file     <- "reference/Homo_sapiens.GRCh38.115.gtf"

out_dir <- "run_bulkRNA/Novel_Microproteins/discovery"
out_prefix   <- file.path(out_dir, "smorf_results")
max_aa_len   <- 100

diamond_evalue    <- "1e-5"
diamond_max_hits   <- 5

# DIAMOND databases
uniprot_full_db <- "reference/uniprot_full.dmnd"
sprot_db        <- "reference/uniprot_sprot.dmnd"
openprot_db     <- "reference/openprot_human.dmnd"
smprot_db       <- "reference/smprot2_human_ribo.dmnd"

# =========================
# Helper functions
# =========================

parse_transdecoder_header <- function(hdr) {
  # Example:
  # ENST00000832824.p1 type:complete gc:universal ENST00000832824:764-345(-)

  first_token <- sub(" .*", "", hdr)

  tx_orf <- first_token
  tx_id <- sub("\\.p[0-9]+$", "", tx_orf)
  orf_id <- tx_orf

  coord_str <- NA_character_
  if (grepl(":[0-9]+-[0-9]+\\([+-]\\)$", hdr)) {
    coord_str <- sub(".* ([^ ]+:[0-9]+-[0-9]+\\([+-]\\))$", "\\1", hdr)
  }

  transcript_id_from_coord <- NA_character_
  tx_start <- NA_integer_
  tx_end <- NA_integer_
  strand <- NA_character_

  if (!is.na(coord_str)) {
    m <- str_match(coord_str, "^(.+):([0-9]+)-([0-9]+)\\(([+-])\\)$")
    transcript_id_from_coord <- m[, 2]
    tx_start <- suppressWarnings(as.integer(m[, 3]))
    tx_end <- suppressWarnings(as.integer(m[, 4]))
    strand <- m[, 5]
  }

  data.table(
    header = hdr,
    orf_id = orf_id,
    tx_id = tx_id,
    tx_id_from_coord = transcript_id_from_coord,
    tx_start = tx_start,
    tx_end = tx_end,
    strand = strand
  )
}

clean_interval <- function(a, b) {
  c(min(a, b), max(a, b))
}

map_tx_interval_to_genome <- function(tx_start, tx_end, exons_df) {
  if (nrow(exons_df) == 0 || any(is.na(c(tx_start, tx_end)))) {
    return(list(
      chrom = NA_character_,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = NA_character_
    ))
  }

  exons_df <- as.data.table(exons_df)

  strand <- unique(exons_df$strand)
  chrom <- unique(as.character(exons_df$seqnames))
  strand <- strand[1]
  chrom <- chrom[1]

  if (strand == "+") {
    setorder(exons_df, start, end)
  } else {
    setorder(exons_df, -start, -end)
  }

  exons_df[, exon_len := end - start + 1L]
  exons_df[, tx_exon_start := cumsum(c(0L, head(exon_len, -1L))) + 1L]
  exons_df[, tx_exon_end := cumsum(exon_len)]

  map_one_pos <- function(tx_pos) {
    hit <- exons_df[tx_exon_start <= tx_pos & tx_exon_end >= tx_pos][1]
    if (nrow(hit) == 0) return(NA_integer_)

    offset <- tx_pos - hit$tx_exon_start
    if (strand == "+") {
      hit$start + offset
    } else {
      hit$end - offset
    }
  }

  g1 <- map_one_pos(tx_start)
  g2 <- map_one_pos(tx_end)

  if (is.na(g1) || is.na(g2)) {
    return(list(
      chrom = chrom,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = strand
    ))
  }

  list(
    chrom = chrom,
    genomic_start = min(g1, g2),
    genomic_end = max(g1, g2),
    strand = strand
  )
}

run_diamond <- function(query_fa, db, out_tsv, evalue = "1e-5", max_hits = 5) {
  diamond_cmd <- c(
    "blastp",
    "-q", query_fa,
    "-d", db,
    "-o", out_tsv,
    "-e", evalue,
    "-k", as.character(max_hits),
    "--quiet",
    "--outfmt", "6",
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )

  message("Running DIAMOND against: ", db)
  status <- system2("diamond", diamond_cmd)
  if (!identical(status, 0L)) {
    warning("DIAMOND returned non-zero exit status for: ", db)
  }
}

read_best_hits <- function(tsv_file, db_name) {
  if (!file.exists(tsv_file) || file.info(tsv_file)$size == 0) return(NULL)

  hits <- fread(tsv_file, header = FALSE)
  setnames(hits, c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  ))

  hits[, `:=`(
    pident = as.numeric(pident),
    length = as.numeric(length),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    qlen = as.numeric(qlen),
    slen = as.numeric(slen)
  )]

  hits[, `:=`(
    qcov = length / qlen,
    scov = length / slen,
    db = db_name
  )]

  # Best hit = lowest evalue, then highest bitscore, then highest identity, then highest qcov
  setorder(hits, qseqid, evalue, -bitscore, -pident, -qcov)
  hits[, .SD[1], by = qseqid]
}

# =========================
# Read and filter peptides
# =========================

pep <- readAAStringSet(pep_fasta)
hdrs <- names(pep)
aa_len <- width(pep)

meta <- rbindlist(lapply(hdrs, parse_transdecoder_header), fill = TRUE)
meta[, aa_len := aa_len]

# Keep only ORFs below the threshold
keep_idx <- aa_len <= max_aa_len
pep_smorf <- pep[keep_idx]
meta_smorf <- meta[keep_idx]

# Write filtered peptide FASTA
smorf_pep_fa <- paste0(out_prefix, ".smorfs_lt", max_aa_len, "aa.pep.fa")
names(pep_smorf) <- meta_smorf$orf_id
writeXStringSet(pep_smorf, smorf_pep_fa)

# =========================
# DIAMOND searches
# =========================

openprot_tsv <- paste0(out_prefix, ".openprot.tsv")
smprot_tsv   <- paste0(out_prefix, ".smprot.tsv")
sprot_tsv    <- paste0(out_prefix, ".uniprot_sprot.tsv")
uniprot_tsv  <- paste0(out_prefix, ".uniprot_full.tsv")

run_diamond(smorf_pep_fa, openprot_db, openprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, smprot_db,   smprot_tsv,   evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, sprot_db,    sprot_tsv,    evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, uniprot_full_db, uniprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)

openprot_best <- read_best_hits(openprot_tsv, "openprot")
smprot_best   <- read_best_hits(smprot_tsv,   "smprot")
sprot_best    <- read_best_hits(sprot_tsv,    "sprot")
uniprot_best  <- read_best_hits(uniprot_tsv,  "uniprot_full")

# =========================
# GTF processing
# =========================

tx_anno <- NULL
tx_exons <- NULL

if (!is.null(gtf_file) && file.exists(gtf_file)) {
  message("Reading GTF...")
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  rm(gtf)
  gc()

  if (!("transcript_id" %in% names(gtf_df))) {
    stop("GTF does not contain a transcript_id column in metadata.")
  }

  anno_cols <- intersect(
    c("transcript_id", "gene_id", "gene_name", "gene_type", "transcript_type", "biotype"),
    names(gtf_df)
  )

  if (length(anno_cols) > 0) {
    tx_anno <- unique(gtf_df[, anno_cols, drop = FALSE])
  }

  exons_only <- gtf_df[gtf_df$type == "exon", , drop = FALSE]

  if (nrow(exons_only) > 0) {
    exons_only <- exons_only[, intersect(
      c("seqnames", "start", "end", "strand", "transcript_id"),
      names(exons_only)
    ), drop = FALSE]

    tx_exons <- as.data.table(exons_only)
    setkey(tx_exons, transcript_id)
  }
}

# =========================
# Build annotation table
# =========================

smorf_dt <- copy(meta_smorf)

# Use transcript ID from coordinate field if available
smorf_dt[, transcript_id := fifelse(!is.na(tx_id_from_coord), tx_id_from_coord, tx_id)]

# Initialize columns
smorf_dt[, `:=`(
  openprot_hit = NA_character_,
  openprot_pident = NA_real_,
  openprot_qcov = NA_real_,
  openprot_evalue = NA_real_,
  openprot_bitscore = NA_real_,

  smprot_hit = NA_character_,
  smprot_pident = NA_real_,
  smprot_qcov = NA_real_,
  smprot_evalue = NA_real_,
  smprot_bitscore = NA_real_,

  sprot_hit = NA_character_,
  sprot_pident = NA_real_,
  sprot_qcov = NA_real_,
  sprot_evalue = NA_real_,
  sprot_bitscore = NA_real_,

  uniprot_hit = NA_character_,
  uniprot_pident = NA_real_,
  uniprot_qcov = NA_real_,
  uniprot_evalue = NA_real_,
  uniprot_bitscore = NA_real_,

  genomic_chrom = NA_character_,
  genomic_start = NA_integer_,
  genomic_end = NA_integer_,
  genomic_strand = NA_character_
)]

# Attach DIAMOND hits by qseqid -> orf_id
if (!is.null(openprot_best) && nrow(openprot_best) > 0) {
  smorf_dt[openprot_best, on = c(orf_id = "qseqid"), `:=`(
    openprot_hit = i.sseqid,
    openprot_pident = i.pident,
    openprot_qcov = i.qcov,
    openprot_evalue = i.evalue,
    openprot_bitscore = i.bitscore
  )]
}

if (!is.null(smprot_best) && nrow(smprot_best) > 0) {
  smorf_dt[smprot_best, on = c(orf_id = "qseqid"), `:=`(
    smprot_hit = i.sseqid,
    smprot_pident = i.pident,
    smprot_qcov = i.qcov,
    smprot_evalue = i.evalue,
    smprot_bitscore = i.bitscore
  )]
}

if (!is.null(sprot_best) && nrow(sprot_best) > 0) {
  smorf_dt[sprot_best, on = c(orf_id = "qseqid"), `:=`(
    sprot_hit = i.sseqid,
    sprot_pident = i.pident,
    sprot_qcov = i.qcov,
    sprot_evalue = i.evalue,
    sprot_bitscore = i.bitscore
  )]
}

if (!is.null(uniprot_best) && nrow(uniprot_best) > 0) {
  smorf_dt[uniprot_best, on = c(orf_id = "qseqid"), `:=`(
    uniprot_hit = i.sseqid,
    uniprot_pident = i.pident,
    uniprot_qcov = i.qcov,
    uniprot_evalue = i.evalue,
    uniprot_bitscore = i.bitscore
  )]
}

# =========================
# Optional genomic mapping
# =========================

if (!is.null(tx_exons)) {
  map_res <- lapply(seq_len(nrow(smorf_dt)), function(i) {
    tx <- smorf_dt$transcript_id[i]
    if (is.na(tx) || !(tx %in% tx_exons$transcript_id)) {
      return(list(
        chrom = NA_character_,
        genomic_start = NA_integer_,
        genomic_end = NA_integer_,
        strand = NA_character_
      ))
    }

    tx_df <- tx_exons[J(tx), nomatch = 0L]
    tx_df <- as.data.frame(tx_df)
    tx_df <- tx_df[, intersect(c("seqnames", "start", "end", "strand"), names(tx_df)), drop = FALSE]

    interval <- clean_interval(smorf_dt$tx_start[i], smorf_dt$tx_end[i])
    map_tx_interval_to_genome(interval[1], interval[2], tx_df)
  })

smorf_dt[, genomic_chrom := vapply(map_res, function(x) as.character(x$chrom), character(1))]
smorf_dt[, genomic_start := vapply(map_res, function(x) as.integer(x$genomic_start), integer(1))]
smorf_dt[, genomic_end := vapply(map_res, function(x) as.integer(x$genomic_end), integer(1))]
smorf_dt[, genomic_strand := vapply(map_res, function(x) as.character(x$strand), character(1))]
}

# =========================
# Classification logic
# =========================

smorf_dt[, `:=`(
  hit_in_openprot = !is.na(openprot_hit),
  hit_in_smprot = !is.na(smprot_hit),

  strong_hit_in_sprot = !is.na(sprot_hit) &
    !is.na(sprot_evalue) & sprot_evalue <= 1e-20 &
    !is.na(sprot_pident) & sprot_pident >= 90 &
    !is.na(sprot_qcov) & sprot_qcov >= 0.80,

  weak_hit_in_uniprot = !is.na(uniprot_hit) &
    !is.na(uniprot_evalue) & uniprot_evalue <= 1e-5 &
    !is.na(uniprot_pident) & uniprot_pident >= 30 &
    !is.na(uniprot_qcov) & uniprot_qcov >= 0.50
)]

smorf_dt[, class := fifelse(
  hit_in_openprot | hit_in_smprot,
  "known_microprotein",
  fifelse(
    strong_hit_in_sprot,
    "known_canonical",
    fifelse(
      weak_hit_in_uniprot,
      "homologous",
      "novel_candidate"
    )
  )
)]

# =========================
# Final cleanup / output
# =========================

# Optional: order columns a little more nicely
setcolorder(smorf_dt, c(
  "header", "orf_id", "tx_id", "tx_id_from_coord", "transcript_id",
  "tx_start", "tx_end", "strand", "aa_len",
  "genomic_chrom", "genomic_start", "genomic_end", "genomic_strand",
  "openprot_hit", "openprot_pident", "openprot_qcov", "openprot_evalue", "openprot_bitscore",
  "smprot_hit", "smprot_pident", "smprot_qcov", "smprot_evalue", "smprot_bitscore",
  "sprot_hit", "sprot_pident", "sprot_qcov", "sprot_evalue", "sprot_bitscore",
  "uniprot_hit", "uniprot_pident", "uniprot_qcov", "uniprot_evalue", "uniprot_bitscore",
  "hit_in_openprot", "hit_in_smprot", "strong_hit_in_sprot", "weak_hit_in_uniprot",
  "class"
))

fwrite(
  smorf_dt,
  paste0(out_prefix, ".smorf_classification.tsv"),
  sep = "\t",
  na = "NA"
)

message("Done. Wrote: ", paste0(out_prefix, ".smorf_classification.tsv"))
```

### Run Novel Microproteins Discovery Configuration File

```bash
cd /data/mckeeka/bulkRNA_sarcoma/

sbatch --job-name=smorf --cpus-per-task=4 --mem=64G --time=8:00:00 --output=smorf_%j.out --error=smorf_%j.err --wrap="set -x; echo START $(date); source ~/.bashrc; conda activate smorf; which Rscript; echo ENV_OK; Rscript Microproteins_Novel_Discovery.r; echo DONE $(date)"

```
