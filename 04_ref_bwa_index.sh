#!/bin/bash
# BWA Index using Docker container
# Output will be in the same directory as the reference genome

# Exit on any error
set -e

# Reference genome path
REF_GENOME="$PWD/01_Data/02_Ref_genome/hg38.fa"

# Directory containing the reference genome
REF_DIR=$(dirname "$REF_GENOME")
REF_BASE=$(basename "$REF_GENOME")

# Docker container with BWA installed
CONTAINER="staphb/bwa:latest"

# Run BWA index inside the container
docker run --rm \
    -v "$REF_DIR":/ref \
    $CONTAINER \
    bash -c "bwa index -p /ref/${REF_BASE} /ref/${REF_BASE}"

echo "BWA index completed. Index files are in $REF_DIR"