#!/usr/bin/env bash
set -euo pipefail
# Discovery additive GWAS. Input is the processed genotype and normalized phenotype set.
CHR=${1:-1}
mapfile -t traits < <(tail -n +2 data/trait_manifest.tsv | cut -f1)
mkdir -p results/gwas
plink2 --bfile "data/genotypes/chr${CHR}" \
    --exclude "data/variant_exclusions/chr${CHR}.txt" --keep data/discovery_keep.tsv \
    --rm-dup force-first --pheno data/gwas_phenotypes.tsv --pheno-name "${traits[@]}" \
    --covar data/gwas_covariates.tsv --maf 0.01 --geno 0.05 --hwe 1e-6 0 \
    --vif 1000 --glm hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,beta,se,tz,p \
    --variance-standardize --no-input-missing-phenotype --seed 9912 \
    --out "results/gwas/chr${CHR}"
