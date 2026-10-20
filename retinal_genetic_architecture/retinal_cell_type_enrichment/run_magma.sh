#!/usr/bin/env bash
set -euo pipefail
# MAGMA v1.10, NCBI37.3 genes, 1000 Genomes EUR; main analysis uses +/-10 kb.
mkdir -p results/scdrs/magma/genes
magma --annotate window=10,10 --snp-loc data/reference/g1000_eur.bim \
    --gene-loc data/reference/NCBI37.3.gene.loc --out results/scdrs/magma/annotation
while read -r trait; do
    magma --bfile data/reference/g1000_eur \
        --pval "results/scdrs/gwas_pval/$trait.pval" use=SNP,P ncol=N \
        --gene-annot results/scdrs/magma/annotation.genes.annot \
        --out "results/scdrs/magma/genes/$trait"
done < <(tail -n +2 data/trait_manifest.tsv | cut -f1)
