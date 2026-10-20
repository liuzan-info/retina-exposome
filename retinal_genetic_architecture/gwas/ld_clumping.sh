#!/usr/bin/env bash
set -euo pipefail
# Construct the discovery LD reference with the original GWAS variant filters.
mkdir -p results/gwas_clump_reference results/gwas_clumps
for chrom in {1..22}; do
    plink2 --bfile "data/genotypes/chr${chrom}" \
        --exclude "data/variant_exclusions/chr${chrom}.txt" --keep data/discovery_keep.tsv \
        --rm-dup force-first --maf 0.01 --geno 0.05 --hwe 1e-6 0 --make-bed \
        --out "results/gwas_clump_reference/chr${chrom}"
done
while read -r trait; do
    mkdir -p "results/gwas_clumps/$trait"
    for chrom in {1..22}; do
        plink2 --bfile "results/gwas_clump_reference/chr${chrom}" \
            --clump "results/gwas/chr${chrom}.${trait}.glm.linear" \
            --clump-p-field P --clump-id-field ID --clump-p1 5e-8 --clump-p2 0.05 \
            --clump-r2 0.1 --clump-kb 1000 --out "results/gwas_clumps/$trait/chr${chrom}"
    done
done < <(tail -n +2 data/trait_manifest.tsv | cut -f1)
