# WGS Variant Calling Pipeline (Nextflow + GATK)

[![Nextflow](https://img.shields.io/badge/Nextflow-v23.04.0-orange)](https://www.nextflow.io/)  
[![GATK](https://img.shields.io/badge/GATK-4.3-blue)](https://gatk.broadinstitute.org/)  
[![Docker](https://img.shields.io/badge/Docker-Enabled-green)](https://www.docker.com/)

- A **reproducible Nextflow pipeline** for calling and annotating SNPs and INDELs from **whole genome sequencing (WGS)** data following **GATK Best Practices**.  
- This pipeline is **Docker-based**, tested on **AWS EC2**, and designed for **reproducibility, scalability, and performance**.

---

## Pipeline Overview

This workflow performs:

1. **Raw read QC**: FastQC + MultiQC  
2. **Read trimming**: Trimmomatic  
3. **Alignment**: BWA-MEM  
4. **Sorting and duplicate marking**: Samtools + Picard  
5. **Base quality score recalibration (BQSR)**  
6. **Variant calling**: GATK HaplotypeCaller  
7. **Variant splitting**: SNPs / INDELs  
8. **Hard filtering of variants**  
9. **Genotype-level filtering**  
10. **Functional annotation**: Funcotator  
11. **Conversion to tabular format**  
12. **Curated result generation**

---

## 📥 Input Requirements

### Raw reads
- Location: `01_Data/01_Raw_Reads/`
- Naming: `*_1.fastq.gz` and `*_2.fastq.gz`

### Reference genome
- Location: `01_Data/02_Ref_genome/`
- Must contain:
  - `hg38.fa`
  - `.fai` index
  - `.dict` sequence dictionary
  - dbSNP VCF

### Funcotator data
- Location: `01_Data/03_Funcotator/`
- Download from:
  https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator

## 📂 Directory Structure

<details>
<summary>Click to see dir structure</summary>

```text
variant_calling_nf/
├── 01_Data
│   ├── 01_Raw_Reads
│   │   ├── SRR062634_1.fastq.gz
│   │   ├── SRR062634_2.fastq.gz
│   │   ├── SRR062644_1.fastq.gz
│   │   └── SRR062644_2.fastq.gz
│   ├── 02_Ref_genome
│   │   ├── Homo_sapiens_assembly38.dbsnp138.vcf
│   │   ├── Homo_sapiens_assembly38.dbsnp138.vcf.idx
│   │   ├── hg38.dict
│   │   ├── hg38.fa
│   │   ├── hg38.fa.amb
│   │   ├── hg38.fa.ann
│   │   ├── hg38.fa.bwt
│   │   ├── hg38.fa.fai
│   │   ├── hg38.fa.pac
│   │   └── hg38.fa.sa
│   └── 03_Funcotator
│      └── funcotator_dataSources.v1.8.hg38.20230908g
├── 01_docker.sh
├── 02_Results
│   ├── 01_Reads_QC
│   │   ├── 01_Raw_fastqc
│   │   │   ├── SRR062634_1_fastqc.html
│   │   │   ├── SRR062634_1_fastqc.zip
│   │   │   ├── SRR062634_2_fastqc.html
│   │   │   ├── SRR062634_2_fastqc.zip
│   │   │   ├── SRR062644_1_fastqc.html
│   │   │   ├── SRR062644_1_fastqc.zip
│   │   │   ├── SRR062644_2_fastqc.html
│   │   │   └── SRR062644_2_fastqc.zip
│   │   ├── 02_MultiQC_Raw_reads
│   │   │   └── multiqc_raw_reads_report.html
│   │   ├── 03_Trimmed_fastqc
│   │   │   ├── SRR062634_trimmed_1_fastqc.html
│   │   │   ├── SRR062634_trimmed_1_fastqc.zip
│   │   │   ├── SRR062634_trimmed_2_fastqc.html
│   │   │   ├── SRR062634_trimmed_2_fastqc.zip
│   │   │   ├── SRR062644_trimmed_1_fastqc.html
│   │   │   ├── SRR062644_trimmed_1_fastqc.zip
│   │   │   ├── SRR062644_trimmed_2_fastqc.html
│   │   │   └── SRR062644_trimmed_2_fastqc.zip
│   │   └── 04_MultiQC_Trimmed_reads
│   │       └── multiqc_trimmed_reads_report.html
│   ├── 02_Trimmed_reads
│   │   ├── SRR062634_trimmed_1.fastq.gz
│   │   ├── SRR062634_trimmed_2.fastq.gz
│   │   ├── SRR062644_trimmed_1.fastq.gz
│   │   └── SRR062644_trimmed_2.fastq.gz
│   ├── 03_Alignment
│   │   ├── 01_BWA_MEM
│   │   │   ├── SRR062634.bam
│   │   │   └── SRR062644.bam
│   │   ├── 02_Sorted
│   │   │   ├── SRR062634.sorted.bam
│   │   │   └── SRR062644.sorted.bam
│   │   ├── 03_Mark_Dup
│   │   │   ├── SRR062634.md.bam
│   │   │   ├── SRR062634.metrics.txt
│   │   │   ├── SRR062644.md.bam
│   │   │   └── SRR062644.metrics.txt
│   │   └── 04_BQSR
│   │       ├── SRR062634.bqsr.bai
│   │       ├── SRR062634.bqsr.bam
│   │       ├── SRR062644.bqsr.bai
│   │       ├── SRR062644.bqsr.bam
│   │       └── recal.table
│   └── 04_Variants
│       ├── 01_Raw
│       │   ├── SRR062634_raw.vcf
│       │   └── SRR062644_raw.vcf
│       ├── 02_Split
│       │   ├── SRR062634_raw_indels.vcf
│       │   ├── SRR062634_raw_snps.vcf
│       │   ├── SRR062644_raw_indels.vcf
│       │   └── SRR062644_raw_snps.vcf
│       ├── 03_Filtered
│       │   ├── SRR062634_filtered_indels.vcf
│       │   ├── SRR062634_filtered_snps.vcf
│       │   ├── SRR062644_filtered_indels.vcf
│       │   └── SRR062644_filtered_snps.vcf
│       ├── 04_Pass
│       │   ├── SRR062634_indels_pass.vcf
│       │   ├── SRR062634_snps_pass.vcf
│       │   ├── SRR062644_indels_pass.vcf
│       │   └── SRR062644_snps_pass.vcf
│       ├── 05_GT_Filtered
│       │   ├── SRR062634_indels_gt_clean.vcf
│       │   ├── SRR062634_snps_gt_clean.vcf
│       │   ├── SRR062644_indels_gt_clean.vcf
│       │   └── SRR062644_snps_gt_clean.vcf
│       ├── 06_Funcotator
│       │   ├── SRR062634_indels_funcotated.vcf
│       │   ├── SRR062634_snps_funcotated.vcf
│       │   ├── SRR062644_indels_funcotated.vcf
│       │   └── SRR062644_snps_funcotated.vcf
│       ├── 07_Tables
│       │   ├── SRR062634_indels.table
│       │   ├── SRR062634_snps.table
│       │   ├── SRR062644_indels.table
│       │   └── SRR062644_snps.table
│       └── 08_Curated_Results
│           ├── SRR062634_curated_indels.txt
│           ├── SRR062634_curated_snps.txt
│           ├── SRR062644_curated_indels.txt
│           └── SRR062644_curated_snps.txt
├── 02_get_demo_data.sh
├── 03_prepare_reference.sh
├── 04_ref_bwa_index.sh
├── 05_Install_JAVA_and_NF.sh
├── 06_Curated_Results.sh
├── nextflow.config
├── variant_calling.nf
├── workflow.html
└── README.md
```
</details>

---

## 🛠 Installation

To install nad run the nextflow pipeline, follow these steps:

1. Clone this repository:

```bash
git clone https://github.com/sgr308/variant_calling_nf.git
```

2. Navigate to the pipeline directory:

```bash
cd variant_calling_nf
```

3. First make all bash scripts executable:
```bash
chmod a+x *.sh
```

4. Run all bash scripts one by one:
```bash
./01_docker.sh # This script will install Docker.
./02_get_demo_data.sh # This script will download one sample of WGS data from 10000 genome project.
./03_prepare_reference.sh # This script will download, index the human reference genonme hg38.fa, creates sequence dictionary and download  known sites (dbSNP).
./04_ref_bwa_index.sh # This script will use BWA to index reference genome.
./05_Install_JAVA_and_NF.sh # This script will install JAVA and nextflow.
```

## Usage
Run the nextflow workflow:
```bash
nextflow run variant_calling.nf -profile docker -with-dag workflow.html
```

## Run Curated Results bash script:
once the nextflow pipeline finishes, run the following script :
```bash
./06_Curated_Results.sh # run this script to get final results
```

## Workflow:

<img src="https://github.com/sgr308/variant_calling_nf/blob/main/nf_variant_calling.jpeg?raw=true"/>
