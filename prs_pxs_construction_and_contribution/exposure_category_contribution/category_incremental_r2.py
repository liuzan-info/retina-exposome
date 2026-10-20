#!/usr/bin/env python3
# Single-category PXS incremental R2 on a common sample; original QR/FWL calculation.
import csv
from pathlib import Path
import numpy as np
import pandas as pd
TOLERANCE = 1e-10
NUMERIC_ZERO_TOLERANCE = 1e-12
MISSING = {"", "NA", "NaN", "nan", "N/A"}

def read_table(path, delimiter):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        columns = {name: [] for name in header}
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{path}: expected {len(header)} fields, found {len(row)}")
            for name, value in zip(header, row):
                columns[name].append(value)
    return header, columns

def numeric(values, name):
    result = np.empty(len(values), dtype=float)
    for index, value in enumerate(values):
        result[index] = np.nan if value in MISSING else float(value)
    if np.isinf(result).any():
        raise ValueError(f"{name}: infinite numeric value")
    return result

def one_hot_sex(values):
    levels = sorted(set(values))
    if not levels:
        raise ValueError("Sex has no observed level")
    return np.column_stack(
        [np.asarray([float(value == level) for value in values]) for level in levels[1:]]
    ) if len(levels) > 1 else np.empty((len(values), 0), dtype=float)

def model_for_trait(trait, feature, covariate, category_scores, available_categories):
    y_all = numeric(feature[trait], f"feature:{trait}")
    age_all = numeric(covariate["Age"], "covariate:Age")
    avmse_all = numeric(covariate["avMSE"], "covariate:avMSE")
    sex_all = np.asarray(covariate["Sex"], dtype=object)
    score_all = {
        category: numeric(category_scores[category][trait], f"{category}:{trait}")
        for category in available_categories
    }

    complete = np.isfinite(y_all) & np.isfinite(age_all) & np.isfinite(avmse_all)
    complete &= np.asarray([value not in MISSING for value in sex_all])
    for score in score_all.values():
        complete &= np.isfinite(score)
    if not complete.any():
        raise ValueError(f"{trait}: no common complete cases")

    y = y_all[complete]
    age = age_all[complete]
    avmse = avmse_all[complete]
    sex_values = sex_all[complete].tolist()
    sex_matrix = one_hot_sex(sex_values)
    baseline = np.column_stack([np.ones(len(y)), age, sex_matrix, avmse])
    baseline_rank = int(np.linalg.matrix_rank(baseline))
    if baseline_rank != baseline.shape[1]:
        raise ValueError(f"{trait}: baseline design is rank deficient")

    q, _ = np.linalg.qr(baseline, mode="reduced")
    y_residual = y - q @ (q.T @ y)
    tss = float(np.sum((y - np.mean(y)) ** 2))
    if not np.isfinite(tss) or tss <= 0:
        raise ValueError(f"{trait}: non-positive outcome TSS")
    rss_base = float(y_residual @ y_residual)
    r2_base = 1.0 - rss_base / tss

    category_rows = []
    raw_scores = []
    for category in available_categories:
        score = score_all[category][complete]
        score_residual = score - q @ (q.T @ score)
        denominator = float(score_residual @ score_residual)
        if denominator <= NUMERIC_ZERO_TOLERANCE:
            raise ValueError(f"{trait} / {category}: category PXS is rank deficient")
        delta = float((score_residual @ y_residual) ** 2 / denominator / tss)
        if delta < -NUMERIC_ZERO_TOLERANCE:
            raise ValueError(f"{trait} / {category}: negative incremental R2 {delta}")
        if abs(delta) < NUMERIC_ZERO_TOLERANCE:
            delta = 0.0
        raw_scores.append(score)
        category_rows.append(
            {
                "EyeFeature": trait,
                "Category": category,
                "n": int(len(y)),
                "r2_base": r2_base,
                "r2_category": r2_base + delta,
                "delta_r2_category": delta,
            }
        )

    delta_sum = sum(row["delta_r2_category"] for row in category_rows)
    for row in category_rows:
        row["category_share"] = (
            row["delta_r2_category"] / delta_sum if delta_sum > 0 else None
        )

    scores = np.column_stack(raw_scores)
    score_residuals = scores - q @ (q.T @ scores)
    score_rank = int(np.linalg.matrix_rank(score_residuals))
    if score_rank != len(available_categories):
        raise ValueError(f"{trait}: joint category PXS design is rank deficient")
    coefficients, _, _, _ = np.linalg.lstsq(score_residuals, y_residual, rcond=None)
    joint_delta = float(y_residual @ (score_residuals @ coefficients) / tss)
    if joint_delta < -NUMERIC_ZERO_TOLERANCE:
        raise ValueError(f"{trait}: negative joint incremental R2 {joint_delta}")
    if abs(joint_delta) < NUMERIC_ZERO_TOLERANCE:
        joint_delta = 0.0

    full_design = np.column_stack([baseline, scores])
    full_rank = int(np.linalg.matrix_rank(full_design))
    if full_rank != baseline.shape[1] + len(available_categories):
        raise ValueError(f"{trait}: direct full model is rank deficient")
    full_coefficients, _, _, _ = np.linalg.lstsq(full_design, y, rcond=None)
    full_residuals = y - full_design @ full_coefficients
    r2_full_direct = 1.0 - float(full_residuals @ full_residuals) / tss
    r2_full = r2_base + joint_delta
    direct_error = r2_full_direct - r2_full
    if abs(direct_error) > TOLERANCE:
        raise ValueError(f"{trait}: FWL/direct full-model error {direct_error}")

    metric = {
        "EyeFeature": trait,
        "n": int(len(y)),
        "n_categories": len(available_categories),
        "available_categories": " | ".join(available_categories),
        "baseline_rank": baseline_rank,
        "full_rank": full_rank,
        "r2_base": r2_base,
        "r2_full_exposome": r2_full,
        "delta_r2_full_exposome": joint_delta,
        "direct_full_model_error": direct_error,
        "interpretation_scope": (
            "Rep1 fitted conditional explained variance; baseline Age,Sex,avMSE; "
            "fixed Step3 category PXS; not causal or cross-validated performance"
        ),
    }
    return category_rows, metric

def build_group_summary(category_rows, included_traits, trait_dictionary, category_dictionary):
    groups = {}
    for trait in included_traits:
        annotation = trait_dictionary[trait]
        group = annotation["group"]
        groups.setdefault(group, []).append(trait)

    lookup = {(row["EyeFeature"], row["Category"]): row for row in category_rows}
    result = []
    for group in sorted(groups, key=lambda value: float(trait_dictionary[groups[value][0]]["group_rank"])):
        group_traits = groups[group]
        raw_means = {}
        coverage = {}
        for category in category_dictionary:
            covered = [trait for trait in group_traits if (trait, category) in lookup]
            shares = [
                lookup[(trait, category)]["category_share"]
                for trait in covered
                if lookup[(trait, category)]["category_share"] is not None
            ]
            deltas = [lookup[(trait, category)]["delta_r2_category"] for trait in covered]
            raw_means[category] = float(np.mean(shares)) if shares else None
            coverage[category] = (covered, shares, deltas)
        normalization = sum(value for value in raw_means.values() if value is not None)
        for category in category_dictionary:
            covered, shares, deltas = coverage[category]
            raw_mean = raw_means[category]
            result.append(
                {
                    "group": group,
                    "group_rank": trait_dictionary[group_traits[0]]["group_rank"],
                    "Category": category,
                    "category_order": category_dictionary[category]["category_order"],
                    "n_group_traits": len(group_traits),
                    "n_covered_traits": len(covered),
                    "n_share_traits": len(shares),
                    "coverage_fraction": len(covered) / len(group_traits),
                    "mean_delta_r2_category_among_covered": float(np.mean(deltas)) if deltas else None,
                    "mean_category_share_among_covered": raw_mean,
                    "renormalized_group_share": (
                        raw_mean / normalization
                        if raw_mean is not None and normalization > 0
                        else None
                    ),
                }
            )
    return result

input_dir = Path("data/scores")
_, feature = read_table(input_dir / "01_rep1_oct_feature.tsv", "\t")
_, covariate = read_table(input_dir / "01_rep1_covariate.tsv", "\t")
assert feature["IID"] == covariate["IID"]
manifest = pd.read_csv(input_dir / "03_group_pxs_z.tsv", sep="\t")
category_scores = {}
available = {}
for row in manifest.itertuples(index=False):
    _, score = read_table(row.path, "\t")
    assert score["IID"] == feature["IID"]
    category_scores[row.group] = score
    available[row.group] = {t for t in score if t != "IID" and any(x not in MISSING for x in score[t])}
traits = sorted((set(feature) - {"IID"}) & set.union(*available.values()))
trait_dictionary = pd.read_csv("data/retinal_dictionary.csv").set_index("EyeFeature").to_dict("index")
category_dictionary = pd.read_csv("data/exposure_categories.csv").set_index("Category").to_dict("index")
category_rows, metrics = [], []
for index, trait in enumerate(traits, 1):
    categories = [c for c in category_scores if trait in available[c]]
    rows, metric = model_for_trait(trait, feature, covariate, category_scores, categories)
    category_rows.extend(rows)
    metrics.append(metric)
    print(f"Category contribution: {index}/{len(traits)} ({100*index/len(traits):.1f}%)", flush=True)
group_rows = build_group_summary(category_rows, traits, trait_dictionary, category_dictionary)
out = Path("results/exposure_category_contribution")
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(category_rows).to_csv(out / "trait_category_incremental_r2.tsv", sep="\t", index=False)
pd.DataFrame(metrics).to_csv(out / "trait_model_metrics.tsv", sep="\t", index=False)
pd.DataFrame(group_rows).to_csv(out / "group_category_incremental_r2.tsv", sep="\t", index=False)
