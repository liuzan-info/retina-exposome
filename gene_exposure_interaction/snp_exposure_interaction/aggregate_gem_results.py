#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("results/vqtl_gem")
RESULT_DIR = ROOT / "results"
ALL_OUT = RESULT_DIR / "gem_interaction_all.tsv"
SIGNIFICANT_OUT = RESULT_DIR / "gem_interaction_significant.tsv"


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(len(ranked), dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def read_gem_result(gem_file: Path) -> pd.DataFrame:
    result = pd.read_csv(gem_file, sep="\t", low_memory=False)
    beta_column = next(column for column in result.columns if column.startswith("Beta_G-"))
    robust_se_column = next(
        column for column in result.columns if column.startswith("robust_SE_Beta_G-")
    )
    result = result.loc[
        :,
        [
            "SNPID",
            "CHR",
            "POS",
            "Non_Effect_Allele",
            "Effect_Allele",
            "N_Samples",
            "AF",
            beta_column,
            robust_se_column,
            "robust_P_Value_Interaction",
        ],
    ].rename(
        columns={
            beta_column: "interaction_beta",
            robust_se_column: "interaction_robust_se",
            "robust_P_Value_Interaction": "interaction_robust_p",
        }
    )
    result.insert(0, "task_id", gem_file.stem)
    result.insert(1, "trait_alias", gem_file.parent.parent.name)
    result.insert(2, "exposure_alias", gem_file.parent.name)
    return result


def main() -> None:
    gem_files = sorted(RESULT_DIR.rglob("*.gem"))
    interactions = pd.concat(
        [read_gem_result(gem_file) for gem_file in gem_files], ignore_index=True
    )

    p_values = interactions["interaction_robust_p"].to_numpy(dtype=float)
    valid = np.isfinite(p_values) & (p_values > 0) & (p_values <= 1)
    valid_count = int(valid.sum())

    interactions["p_fdr_bh"] = np.nan
    interactions["bonferroni_threshold"] = np.nan
    interactions["bonferroni_pass"] = False
    if valid_count:
        valid_p = p_values[valid]
        interactions.loc[valid, "p_fdr_bh"] = bh_adjust(valid_p)
        interactions.loc[valid, "bonferroni_threshold"] = 0.05 / valid_count
        interactions.loc[valid, "bonferroni_pass"] = valid_p <= 0.05 / valid_count

    interactions.to_csv(ALL_OUT, sep="\t", index=False, na_rep="NA")
    significant = interactions.loc[
        valid & ((interactions["p_fdr_bh"] <= 0.05) | interactions["bonferroni_pass"])
    ].copy()
    significant.to_csv(SIGNIFICANT_OUT, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
