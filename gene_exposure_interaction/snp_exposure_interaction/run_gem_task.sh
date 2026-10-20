#!/usr/bin/env bash

set -euo pipefail

ROOT="results/vqtl_gem"
GEM_BIN="GEM"
PGEN_PREFIX="$ROOT/genotypes/discovery_lead_vqtl"
MASTER_PHENO="$ROOT/inputs/discovery_gem_master.csv"
COVAR_NAMES_FILE="$ROOT/inputs/gem_covar_names.txt"
ANALYSIS_MANIFEST="$ROOT/manifests/analysis_manifest_all.tsv"
RESULT_DIR="$ROOT/results"

TASK_ID="$1"
THREADS="${2:-5}"
record=$(awk -F '\t' -v task="$TASK_ID" 'NR > 1 && $1 == task {print $3 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $11; exit}' "$ANALYSIS_MANIFEST")
IFS=$'\t' read -r trait_alias trait_column exposure_alias exposure_column exposure_type snp_file <<< "$record"

output_dir="$RESULT_DIR/$trait_alias/$exposure_alias"
mkdir -p "$output_dir"
covariates=$(paste -sd ' ' "$COVAR_NAMES_FILE")

gem_args=(
    --pfile "$PGEN_PREFIX"
    --pheno-file "$MASTER_PHENO"
    --sampleid-name IID
    --pheno-name "$trait_column"
    --exposure-names "$exposure_column"
    --covar-names $covariates
    --include-snp-file "$snp_file"
    --delim ,
    --maf 0.05
    --missing-value NA
    --robust 1
    --cat-threshold 2
)
if [[ "$exposure_type" == "Binary" || "$exposure_type" == "CatUnord" ]]; then
    gem_args+=(--categorical-names "$exposure_column")
fi
gem_args+=(
    --threads "$THREADS"
    --stream-snps 100
    --output-style meta
    --out "$output_dir/$TASK_ID.gem"
)

OMP_NUM_THREADS="$THREADS" "$GEM_BIN" "${gem_args[@]}"
