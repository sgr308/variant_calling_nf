#!/usr/bin/env bash

set -e  # Exit on error

echo "Checking if Docker is installed..."

if command -v docker >/dev/null 2>&1; then
    echo "Docker is already installed."
else
    echo "Docker not found. Installing Docker..."

    # Update package list
    sudo apt update

    # Install Docker
    sudo apt install -y docker.io
fi

echo "Starting Docker service..."
sudo systemctl start docker
sudo systemctl enable docker

echo "Docker version:"
docker --version
