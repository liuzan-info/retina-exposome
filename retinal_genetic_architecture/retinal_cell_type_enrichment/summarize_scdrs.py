#!/usr/bin/env python3
# Final main-cell-type statistics; no plotting or donor bootstrap.
from pathlib import Path
import numpy as np
import pandas as pd
CELL_TYPES = ["rod", "cone", "bipolar", "amacrine", "horizontal", "RGC", "Muller_glia", "astrocyte", "microglia", "RPE"]
def bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(adjusted, 1.0)
    return result

def read_group_results(root, annotation, seeds):
    rows = []
    for seed in seeds:
        seed_dir = Path(root) / f"seed{seed}"
        for path in sorted(seed_dir.glob(f"*.scdrs_group.{annotation}")):
            trait = path.name.split(".scdrs_group.")[0]
            frame = pd.read_csv(path, sep="\t", index_col=0).reset_index()
            frame = frame.rename(columns={frame.columns[0]: "cell_group"})
            frame["trait"] = trait
            frame["seed"] = seed
            rows.append(frame)
    return pd.concat(rows, ignore_index=True)

def build_final(raw, metadata, final_seeds, min_positive_seeds):
    raw = raw.loc[
        raw["trait"].isin(metadata["trait"])
        & raw["cell_group"].isin(CELL_TYPES)
        & raw["seed"].isin(final_seeds)
    ].copy()
    assert len(raw) == 58 * 10 * len(final_seeds)
    assert not raw.duplicated(["trait", "cell_group", "seed"]).any()
    assert raw["trait"].nunique() == 58
    assert raw["cell_group"].nunique() == 10
    assert set(raw["seed"]) == set(final_seeds)

    raw["prop_fdr_0.1"] = raw["n_fdr_0.1"] / raw["n_cell"]
    result = (
        raw.groupby(["trait", "cell_group"], as_index=False)
        .agg(
            **{
                "prop_fdr_0.1": ("prop_fdr_0.1", "median"),
                "n_seed_available": ("seed", "nunique"),
                "assoc_positive_seed_count": (
                    "assoc_mcz",
                    lambda x: int((x > 0).sum()),
                ),
                "hetero_positive_seed_count": (
                    "hetero_mcz",
                    lambda x: int((x > 0).sum()),
                ),
            }
        )
    )

    seed1 = raw.loc[raw["seed"].eq(1)].merge(
        metadata[["trait", "analysis_tier"]], on="trait", how="left"
    )
    assert seed1["analysis_tier"].value_counts().to_dict() == {
        "exploratory": 300,
        "validation": 280,
    }
    seed1["assoc_fdr"] = np.nan
    seed1["hetero_fdr"] = np.nan
    for _, indices in seed1.groupby("analysis_tier").groups.items():
        seed1.loc[indices, "assoc_fdr"] = bh_adjust(seed1.loc[indices, "assoc_mcp"])
        seed1.loc[indices, "hetero_fdr"] = bh_adjust(
            seed1.loc[indices, "hetero_mcp"]
        )

    result = result.merge(
        seed1[["trait", "cell_group", "assoc_fdr", "hetero_fdr"]],
        on=["trait", "cell_group"],
        how="left",
    )
    result["association_supported"] = (
        result["assoc_fdr"].lt(0.05)
        & result["assoc_positive_seed_count"].ge(min_positive_seeds)
    )
    result["heterogeneity_supported"] = (
        result["association_supported"]
        & result["hetero_fdr"].lt(0.05)
        & result["hetero_positive_seed_count"].ge(min_positive_seeds)
    )
    assert not (
        result["heterogeneity_supported"] & ~result["association_supported"]
    ).any()
    return result

manifest = pd.read_csv("data/trait_manifest.tsv", sep="\t")
raw = read_group_results("results/scdrs/downstream", "major_cell_type", [1, 2, 4, 5])
result = build_final(raw, manifest, [1, 2, 4, 5], 3)
Path("results/scdrs/summary").mkdir(parents=True, exist_ok=True)
raw.to_csv("results/scdrs/summary/major_cell_type_all_seeds.tsv", sep="\t", index=False)
result.to_csv("results/scdrs/summary/major_cell_type_summary.tsv", sep="\t", index=False)
