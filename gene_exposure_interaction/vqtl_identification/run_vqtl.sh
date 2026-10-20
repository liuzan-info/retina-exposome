#!/usr/bin/env bash
set -euo pipefail
mkdir -p results/vqtl/qc_genotype
for chrom in {1..22}; do
    plink2 --bfile "data/genotypes/chr${chrom}" \
        --exclude "data/variant_exclusions/chr${chrom}.txt" --keep data/discovery_keep.tsv \
        --rm-dup force-first --maf 0.05 --mac 50 --geno 0.05 --hwe 1e-6 0 --make-bed \
        --out "results/vqtl/qc_genotype/chr${chrom}"
done
mapfile -t traits < results/vqtl/pheno/phenotype_list.txt
completed=0
total=$((${#traits[@]} * 22))
for trait in "${traits[@]}"; do
    mkdir -p "results/vqtl/scan/$trait" "results/vqtl/clumps/$trait"
    for chrom in {1..22}; do
        osca --vqtl --vqtl-mtd Levene_median --bfile "results/vqtl/qc_genotype/chr${chrom}" \
            --pheno "results/vqtl/pheno/$trait.pheno" --thread-num 1 \
            --out "results/vqtl/scan/$trait/chr${chrom}"
        plink2 --bfile "results/vqtl/qc_genotype/chr${chrom}" \
            --clump "results/vqtl/scan/$trait/chr${chrom}.vqtl" --clump-p-field P --clump-id-field SNP \
            --clump-p1 5e-8 --clump-p2 0.05 --clump-r2 0.1 --clump-kb 1000 \
            --out "results/vqtl/clumps/$trait/chr${chrom}"
        completed=$((completed + 1))
        printf 'vQTL scan/clump: %s/%s\n' "$completed" "$total"
    done
done
