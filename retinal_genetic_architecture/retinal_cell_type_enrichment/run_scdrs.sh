#!/usr/bin/env bash
set -euo pipefail
# Build the 1,000-gene Z-score-weighted sets from the MAGMA gene Z-score matrix.
scdrs munge-gs --out-file results/scdrs/main_10kb.gs \
    --zscore-file results/scdrs/main_10kb.zscore.tsv --weight zscore --n-max 1000
for seed in 1 2 4 5; do
    score_out="results/scdrs/scores/seed$seed"
    downstream_out="results/scdrs/downstream/seed$seed"
    mkdir -p "$score_out" "$downstream_out"
    scdrs compute-score --h5ad-file "data/hrca/hrca_seed${seed}.h5ad" \
        --h5ad-species human --gs-file results/scdrs/main_10kb.gs --gs-species human \
        --out-folder "$score_out" --cov-file "data/hrca/hrca_seed${seed}.cov" \
        --flag-filter-data True --flag-raw-count True --n-ctrl 1000 \
        --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True
    scdrs perform-downstream --h5ad-file "data/hrca/hrca_seed${seed}.h5ad" \
        --score-file "$score_out/@.full_score.gz" --out-folder "$downstream_out" \
        --group-analysis major_cell_type,major_cell_subtype,region \
        --flag-filter-data True --flag-raw-count True
done
