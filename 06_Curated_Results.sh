#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$PWD/02_Results/04_Variants"

FUNC_DIR="${BASE_DIR}/06_Funcotator"
TABLE_DIR="${BASE_DIR}/07_Tables"
OUT_DIR="${BASE_DIR}/08_Curated_Results"

mkdir -p "$OUT_DIR"

# Function to clean table data
clean_table() {
    local table_file="$1"
    # Remove FUNCOTATION header, remove brackets, split glued records, replace pipes with tabs
    sed '/^FUNCOTATION/d' "$table_file" \
      | tr -d '[]' \
      | sed 's/,[[:space:]]*\[/\n/g' \
      | sed 's/|/\t/g'
}

# ----------------------
# Process SNPs
# ----------------------
for vcf in "$FUNC_DIR"/*_snps_funcotated.vcf; do
    id="$(basename "$vcf" _snps_funcotated.vcf)"
    echo "Curating SNPs for $id"

    # Extract Funcotator header cleanly and replace pipes with tabs
    grep "Funcotation fields are:" "$vcf" \
      | sed 's/.*Funcotation fields are: //' \
      | sed 's/\">.*//' \
      | sed 's/|/\t/g' \
      > "${OUT_DIR}/${id}_curated_snps.txt"

    # Append cleaned table data (column 5) with pipes converted to tabs
    cut -f5 "${TABLE_DIR}/${id}_snps.table" \
      | clean_table /dev/stdin \
      >> "${OUT_DIR}/${id}_curated_snps.txt"
done

# ----------------------
# Process INDELs
# ----------------------
for vcf in "$FUNC_DIR"/*_indels_funcotated.vcf; do
    id="$(basename "$vcf" _indels_funcotated.vcf)"
    echo "Curating INDELs for $id"

    # Extract Funcotator header cleanly and replace pipes with tabs
    grep "Funcotation fields are:" "$vcf" \
      | sed 's/.*Funcotation fields are: //' \
      | sed 's/\">.*//' \
      | sed 's/|/\t/g' \
      > "${OUT_DIR}/${id}_curated_indels.txt"

    # Append cleaned table data (column 5) with pipes converted to tabs
    cut -f5 "${TABLE_DIR}/${id}_indels.table" \
      | clean_table /dev/stdin \
      >> "${OUT_DIR}/${id}_curated_indels.txt"
done

echo "All variant curation complete."