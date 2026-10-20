#!/usr/bin/env bash
set -euo pipefail
trait="$1"
n_gwas=27948
seed=1206
chroms="22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1"
mkdir -p "results/prs/weights/$trait"
python software/PRScs/PRScs.py --ref_dir=data/reference/ldblk_1kg_eur \
    --bim_prefix=results/prs/target_allchr --sst_file="results/prs/sumstats/$trait.prscs.txt" \
    --n_gwas="$n_gwas" --chrom="$chroms" --seed="$seed" \
    --out_dir="results/prs/weights/$trait/$trait"
