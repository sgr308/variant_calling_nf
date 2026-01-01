#!/usr/bin/env bash
set -euo pipefail

# Create 01_Data/01_Raw_Reads only if it does not already exist in $PWD
[ -d "01_Data/01_Raw_Reads" ] || mkdir -p "01_Data/01_Raw_Reads"

REFDIR="01_Data/02_Ref_genome"
mkdir -p ${REFDIR}
cd ${REFDIR}

echo "Downloading hg38 reference..."
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -f hg38.fa.gz

echo "Indexing reference with samtools..."
sudo docker run --rm -v $PWD:/data biocontainers/samtools:v1.9-4-deb_cv1 \
    samtools faidx /data/hg38.fa

echo "Creating sequence dictionary..."
sudo docker run --rm -v $PWD:/data broadinstitute/gatk:latest \
    gatk CreateSequenceDictionary \
    -R /data/hg38.fa \
    -O /data/hg38.dict

echo "Downloading known sites (dbSNP)..."
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

echo "Reference preparation complete"

echo "Downloading Funcotator..."
sudo docker run --rm -v $PWD:/data broadinstitute/gatk:latest \
gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download --hg38

echo "Funcotator Downloading complete"
