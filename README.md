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
├── 01_Data/
│   ├── 01_Raw_Reads/          # Raw FASTQ files
│   ├── 02_Ref_genome/         # Reference genome + indexes + known sites
│   │   ├── hg38.fa
│   │   ├── hg38.fa.fai
│   │   ├── hg38.dict
│   │   └── dbsnp.vcf
│   └── 03_Funcotator/         # Funcotator reference data
├── 02_Results/                # Output results (FastQC → Curated tables)
├── 01_docker.sh
├── 02_get_demo_data.sh
├── 03_prepare_reference.sh
├── 04_ref_bwa_index.sh
├── 05_Install_JAVA_and_NF.sh
├── 06_Curated_Results.sh
├── variant_calling.nf         # Main Nextflow workflow
├── nextflow.config            # Nextflow configuration
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
