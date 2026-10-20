#!/usr/bin/env bash
set -euo pipefail
# European LDSC reference and HapMap3 harmonization; only the manuscript h2 branch.
LDSC_DIR="software/ldsc"
mkdir -p results/ldsc/munge results/ldsc/h2
while read -r trait; do
    python "$LDSC_DIR/munge_sumstats.py" --sumstats "data/gwas_summary/${trait}.glm.linear" \
        --snp SNP --N-col OBS_CT --merge-alleles data/reference/w_hm3.snplist \
        --out "results/ldsc/munge/$trait"
    python "$LDSC_DIR/ldsc.py" --h2 "results/ldsc/munge/$trait.sumstats.gz" \
        --ref-ld-chr data/reference/eur_w_ld_chr/ --w-ld-chr data/reference/eur_w_ld_chr/ \
        --out "results/ldsc/h2/$trait"
done < <(tail -n +2 data/trait_manifest.tsv | cut -f1)
extract_h2_summary() {
    local merge_mode="$1"
    local h2_dir="$2"
    local summary_file="$h2_dir/heritability_summary.tsv"

    printf "merge_mode\ttrait\th2_obs\th2_obs_se\tlambda_gc\tmean_chi2\tintercept\tintercept_se\tratio\tratio_se\tlog_file\n" > "$summary_file"

    for log_file in "$h2_dir"/*.log; do
        trait="$(basename "$log_file" .log)"
        awk -v merge_mode="$merge_mode" -v trait="$trait" -v log_file="$log_file" '
            BEGIN {
                h2_obs = "NA"
                h2_obs_se = "NA"
                lambda_gc = "NA"
                mean_chi2 = "NA"
                intercept = "NA"
                intercept_se = "NA"
                ratio = "NA"
                ratio_se = "NA"
            }
            /^Total Observed scale h2:/ {
                value = $0
                sub(/^Total Observed scale h2: /, "", value)
                h2_obs = value
                sub(/ .*/, "", h2_obs)
                h2_obs_se = value
                sub(/^[^(]*\(/, "", h2_obs_se)
                sub(/\).*/, "", h2_obs_se)
            }
            /^Lambda GC:/ {
                lambda_gc = $3
            }
            /^Mean Chi\^2:/ {
                mean_chi2 = $3
            }
            /^Intercept:/ {
                if ($2 == "constrained") {
                    intercept = $0
                    sub(/^Intercept: /, "", intercept)
                    intercept_se = "NA"
                } else {
                    value = $0
                    sub(/^Intercept: /, "", value)
                    intercept = value
                    sub(/ .*/, "", intercept)
                    intercept_se = value
                    sub(/^[^(]*\(/, "", intercept_se)
                    sub(/\).*/, "", intercept_se)
                }
            }
            /^Ratio:/ {
                value = $0
                sub(/^Ratio: /, "", value)
                ratio = value
                sub(/ .*/, "", ratio)
                if (value ~ /\(/) {
                    ratio_se = value
                    sub(/^[^(]*\(/, "", ratio_se)
                    sub(/\).*/, "", ratio_se)
                }
            }
            /^Ratio < 0/ {
                ratio = $0
                ratio_se = "NA"
            }
            END {
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", merge_mode, trait, h2_obs, h2_obs_se, lambda_gc, mean_chi2, intercept, intercept_se, ratio, ratio_se, log_file
            }
        ' "$log_file"
    done >> "$summary_file"
}

extract_h2_summary with_merge_alleles results/ldsc/h2
