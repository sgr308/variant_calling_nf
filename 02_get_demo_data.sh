#!/usr/bin/env bash
set -euo pipefail

mkdir -p 01_Data/01_Raw_Reads
cd 01_Data/01_Raw_Reads

echo "Downloading demo FASTQ files..."

wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

echo "Demo FASTQ download complete"
