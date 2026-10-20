#!/usr/bin/env bash
set -euo pipefail
# Run from repository root. Manifest: trait column; GWAS files: PLINK cols=chrom,pos,ref,alt,a1freq,nobs,beta,se,tz,p.
mkdir -p results/scdrs/gwas_pval
while read -r trait; do
    out="results/scdrs/gwas_pval/$trait.pval"
    printf 'CHR\tBP\tSNP\tP\tN\n' > "$out"
    for chrom in $(seq 1 22); do
        awk 'BEGIN { FS=OFS="\t" } NR > 1 && $3 != "." && $12 != "NA" && $12 != "" { print $1,$2,$3,$12,$8 }' \
            "data/gwas/chr${chrom}.${trait}.glm.linear" >> "$out"
    done
done < <(tail -n +2 data/trait_manifest.tsv | cut -f1)
