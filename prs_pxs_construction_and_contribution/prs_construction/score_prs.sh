#!/usr/bin/env bash
set -euo pipefail
# Usage: bash score_prs.sh rep1  OR  bash score_prs.sh white_non_discovery
# Inputs are preprocessed target genotypes for that population; no discovery samples are scored.
population="$1"
plink2_path="plink2"
target_dir="data/prs_targets/$population"
weight_root="results/prs/weights"
out_dir="results/prs/scores/$population"
per_chr_dir="${out_dir}/per_chr"
per_trait_dir="${out_dir}/per_trait"
matrix_file="${out_dir}/prscs_sum_avg_matrix.tsv"

mapfile -t traits < <(find "$weight_root" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

mkdir -p "$per_chr_dir" "$per_trait_dir"

for trait in "${traits[@]}"; do
    trait_chr_dir="${per_chr_dir}/${trait}"
    mkdir -p "$trait_chr_dir"

    echo "Start PLINK2 scoring ${trait}: $(date)"
    for chr in {1..22}; do
        weight_file="${weight_root}/${trait}/${trait}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

        "$plink2_path" \
            --bfile "${target_dir}/chr${chr}" \
            --score "$weight_file" 2 4 6 cols=fid,nallele,denom,scoreavgs,scoresums \
            --out "${trait_chr_dir}/chr${chr}"
    done

    echo "Summarize chromosome scores ${trait}: $(date)"
    awk '
        BEGIN {
            FS = OFS = "\t"
        }
        FNR == 1 {
            file_idx++
            if (file_idx == 1) {
                for (i = 1; i <= NF; i++) {
                    name = $i
                    sub(/^#/, "", name)
                    col[name] = i
                }
            }
            next
        }
        {
            key = $1 SUBSEP $2
            if (!(key in seen)) {
                seen[key] = 1
                order[++n] = key
                fid[key] = $1
                iid[key] = $2
            }
            sum_score[key] += $(col["SCORE1_SUM"]) + 0
            denom[key] += $(col["DENOM"]) + 0
        }
        END {
            print "FID", "IID", "SUM", "AVG"
            for (i = 1; i <= n; i++) {
                key = order[i]
                avg = sum_score[key] / denom[key]
                printf "%s\t%s\t%.12g\t%.12g\n", fid[key], iid[key], sum_score[key], avg
            }
        }
    ' "${trait_chr_dir}"/chr{1..22}.sscore > "${per_trait_dir}/${trait}.sum_avg.tsv"

    echo "Finished ${trait}: $(date)"
done

trait_names="${traits[*]}"
per_trait_files=()
for trait in "${traits[@]}"; do
    per_trait_files+=("${per_trait_dir}/${trait}.sum_avg.tsv")
done

awk -v trait_names="$trait_names" -v ntraits="${#traits[@]}" '
    BEGIN {
        FS = OFS = "\t"
        split(trait_names, traits, " ")
    }
    FNR == 1 {
        file_idx++
        trait = traits[file_idx]
        next
    }
    {
        key = $1 SUBSEP $2
        if (file_idx == 1) {
            order[++n] = key
            fid[key] = $1
            iid[key] = $2
        }
        sum_score[key, trait] = $3
        avg_score[key, trait] = $4
    }
    END {
        printf "FID\tIID"
        for (i = 1; i <= ntraits; i++) {
            printf "\t%s_SUM\t%s_AVG", traits[i], traits[i]
        }
        printf "\n"

        for (row = 1; row <= n; row++) {
            key = order[row]
            printf "%s\t%s", fid[key], iid[key]
            for (i = 1; i <= ntraits; i++) {
                trait = traits[i]
                printf "\t%s\t%s", sum_score[key, trait], avg_score[key, trait]
            }
            printf "\n"
        }
    }
' "${per_trait_files[@]}" > "$matrix_file"

echo "Final PRS matrix: $matrix_file"
echo "Step3 finished: $(date)"
