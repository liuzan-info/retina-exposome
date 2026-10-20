#!/usr/bin/env bash
set -euo pipefail
# PRS-CS input: SNP A1 A2 BETA SE; target BIM concatenated over 22 autosomes.
target_dir="data/prs_targets/rep1"
summary_dir="data/gwas_summary"
bim_prefix="results/prs/target_allchr"
sumstats_out_dir="results/prs/sumstats"
mkdir -p "$sumstats_out_dir"
: > "${bim_prefix}.bim"
for chr in {1..22}; do
    cat "${target_dir}/chr${chr}.bim" >> "${bim_prefix}.bim"
done

for summary_file in "$summary_dir"/*.glm.linear; do
    summary_name="$(basename "$summary_file" .glm.linear)"
    out_file="${sumstats_out_dir}/${summary_name}.prscs.txt"

    awk 'BEGIN { OFS = "\t" }
         NR == 1 {
             for (i = 1; i <= NF; i++) {
                 col[$i] = i
             }
             print "SNP", "A1", "A2", "BETA", "SE"
             next
         }
         {
             print $(col["SNP"]), $(col["A1"]), $(col["A2"]), $(col["BETA"]), $(col["SE"])
         }' "$summary_file" > "$out_file"
done

