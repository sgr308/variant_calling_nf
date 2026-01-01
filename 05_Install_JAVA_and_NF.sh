#!/usr/bin/env bash
set -euo pipefail

echo "======================================"
echo " Installing Java 17 and Nextflow"
echo "======================================"

# -----------------------------
# Check OS
# -----------------------------
if ! grep -qi "amazon linux" /etc/os-release; then
    echo "WARNING: This script is optimized for Amazon Linux."
fi

# -----------------------------
# Update system
# -----------------------------
sudo yum update -y

# -----------------------------
# Install Java 17 (Amazon Corretto)
# -----------------------------
echo "Installing Java 17..."
sudo yum install -y java-17-amazon-corretto

# -----------------------------
# Verify Java
# -----------------------------
echo "Java version:"
java -version

# -----------------------------
# Install Nextflow
# -----------------------------
echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash

# -----------------------------
# Move Nextflow to PATH
# -----------------------------
sudo mv nextflow /usr/local/bin/

# -----------------------------
# Verify Nextflow
# -----------------------------
echo "Nextflow version:"
nextflow -version

echo "======================================"
echo " Installation complete"
echo "======================================"
