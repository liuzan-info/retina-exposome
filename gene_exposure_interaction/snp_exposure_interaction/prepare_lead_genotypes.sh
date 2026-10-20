#!/usr/bin/env bash

set -euo pipefail

ROOT="results/vqtl_gem"
PLINK_BIN="plink2"
BFILE_PREFIX="results/vqtl/qc_genotype/chr"
KEEP_FILE="data/discovery_keep.tsv"
MANIFEST_DIR="$ROOT/manifests"
GENOTYPE_DIR="$ROOT/genotypes"
PGEN_PREFIX="$GENOTYPE_DIR/discovery_lead_vqtl"
PMERGE_LIST="$GENOTYPE_DIR/pmerge_discovery_lead_vqtl.list"

mkdir -p "$GENOTYPE_DIR"
: > "$PMERGE_LIST"

for chromosome in {1..22}; do
    extract_file="$MANIFEST_DIR/vqtl_union_chr$chromosome.extract"
    output_prefix="$GENOTYPE_DIR/chr$chromosome.discovery_lead_vqtl"

    [[ -s "$extract_file" ]] || continue

    "$PLINK_BIN" \
        --bfile "$BFILE_PREFIX$chromosome" \
        --keep "$KEEP_FILE" \
        --extract "$extract_file" \
        --rm-dup force-first \
        --sort-vars \
        --make-pgen \
        --threads 100 \
        --out "$output_prefix"

    printf '%s\n' "$output_prefix" >> "$PMERGE_LIST"
done

"$PLINK_BIN" \
    --pmerge-list "$PMERGE_LIST" pfile \
    --make-pgen \
    --threads 100 \
    --out "$PGEN_PREFIX"
