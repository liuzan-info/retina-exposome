#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd


ROOT = Path("results/vqtl_gem")
VQTL_CLUMPS = Path("results/vqtl/all_58phenotypes.clumps")
STEP3_EXPOSURE_DICT = Path("data/exposure_dictionary.csv")
KEEP_FILE = Path("data/discovery_keep.tsv")
PHENO_FILE = Path("data/gwas_phenotypes.tsv")
COVAR_FILE = Path("data/gwas_covariates.tsv")
EXPOSURE_FILE = Path("data/discovery_exposures.csv")

MANIFEST_DIR = ROOT / "manifests"
TRAIT_SNPLIST_DIR = MANIFEST_DIR / "trait_snplists"
INPUT_DIR = ROOT / "inputs"
P_THRESHOLD = 5e-8


def trait_sort_key(value: str) -> tuple[int, str]:
    first = value.split("_", 1)[0]
    try:
        return int(first), value
    except ValueError:
        return sys.maxsize, value


def normalize_iid(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)


def main() -> None:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    TRAIT_SNPLIST_DIR.mkdir(parents=True, exist_ok=True)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    clumps = pd.read_csv(
        VQTL_CLUMPS,
        sep="\t",
        dtype={"safe_phenotype": "string", "#CHROM": "string", "ID": "string"},
        low_memory=False,
    )
    leads = clumps.loc[:, ["safe_phenotype", "#CHROM", "POS", "ID", "P"]].copy()
    leads.columns = ["trait", "chromosome", "position", "snp_id", "p_vqtl"]
    leads["p_vqtl"] = pd.to_numeric(leads["p_vqtl"])
    leads["position"] = pd.to_numeric(leads["position"])
    leads = leads.loc[
        leads["trait"].notna()
        & leads["chromosome"].notna()
        & leads["snp_id"].notna()
        & (leads["snp_id"] != ".")
        & leads["p_vqtl"].notna()
        & (leads["p_vqtl"] <= P_THRESHOLD)
    ].copy()
    leads["trait"] = leads["trait"].astype(str)
    leads["chromosome"] = leads["chromosome"].astype(str)
    leads["snp_id"] = leads["snp_id"].astype(str)
    leads = leads.drop_duplicates(["trait", "chromosome", "snp_id"])
    leads = leads.sort_values(
        ["trait", "chromosome", "position", "snp_id"],
        key=lambda column: (
            column.map(trait_sort_key)
            if column.name == "trait"
            else pd.to_numeric(column)
            if column.name in {"chromosome", "position"}
            else column
        ),
        kind="stable",
    ).reset_index(drop=True)

    exposures = pd.read_csv(STEP3_EXPOSURE_DICT, dtype="string", low_memory=False)
    exposures = exposures.loc[
        exposures["expression_exposure"].notna()
        & exposures["association_exposure"].notna()
        & exposures["Type"].notna()
    ].drop_duplicates("expression_exposure", keep="first").reset_index(drop=True)

    traits = sorted(leads["trait"].drop_duplicates().tolist(), key=trait_sort_key)
    trait_alias = {trait: f"P{index:03d}" for index, trait in enumerate(traits, start=1)}
    exposure_names = exposures["expression_exposure"].astype(str).tolist()
    exposure_alias = {
        exposure: f"E{index:03d}" for index, exposure in enumerate(exposure_names, start=1)
    }

    leads.insert(1, "trait_alias", leads["trait"].map(trait_alias))
    leads.to_csv(MANIFEST_DIR / "lead_vqtl_manifest.tsv", sep="\t", index=False, na_rep="NA")

    trait_rows = []
    for trait in traits:
        alias = trait_alias[trait]
        snp_file = TRAIT_SNPLIST_DIR / f"{alias}.snpid"
        snp_ids = leads.loc[leads["trait"] == trait, "snp_id"].drop_duplicates().astype(str)
        pd.DataFrame({"snpid": snp_ids}).to_csv(snp_file, sep="\t", index=False)
        trait_rows.append(
            {
                "trait": trait,
                "trait_alias": alias,
                "lead_snp_count_requested": len(snp_ids),
                "lead_snp_file": str(snp_file),
            }
        )
    trait_manifest = pd.DataFrame(trait_rows)
    trait_manifest.to_csv(MANIFEST_DIR / "trait_manifest.tsv", sep="\t", index=False)

    exposure_manifest = exposures.copy()
    exposure_manifest.insert(
        0,
        "exposure_alias",
        exposure_manifest["expression_exposure"].astype(str).map(exposure_alias),
    )
    exposure_manifest.to_csv(MANIFEST_DIR / "exposure_manifest.tsv", sep="\t", index=False)

    union_ids = sorted(leads["snp_id"].drop_duplicates().astype(str).tolist())
    (MANIFEST_DIR / "vqtl_union.extract").write_text(
        "\n".join(union_ids) + "\n", encoding="utf-8"
    )
    for chromosome in range(1, 23):
        (MANIFEST_DIR / f"vqtl_union_chr{chromosome}.extract").write_text(
            "", encoding="utf-8"
        )
    for chromosome, subset in leads.groupby("chromosome", sort=False):
        ids = sorted(subset["snp_id"].drop_duplicates().astype(str).tolist())
        (MANIFEST_DIR / f"vqtl_union_chr{chromosome}.extract").write_text(
            "\n".join(ids) + "\n", encoding="utf-8"
        )

    task_rows = []
    for trait_row in trait_rows:
        for exposure in exposure_manifest.itertuples(index=False):
            task_rows.append(
                {
                    "task_id": f"{trait_row['trait_alias']}_{exposure.exposure_alias}",
                    "trait": trait_row["trait"],
                    "trait_alias": trait_row["trait_alias"],
                    "trait_column": trait_row["trait_alias"],
                    "expression_exposure": exposure.expression_exposure,
                    "association_exposure": exposure.association_exposure,
                    "exposure_alias": exposure.exposure_alias,
                    "exposure_column": exposure.exposure_alias,
                    "exposure_type": exposure.Type,
                    "semantic_numeric_opposite": exposure.semantic_numeric_opposite,
                    "lead_snp_file": trait_row["lead_snp_file"],
                }
            )
    pd.DataFrame(task_rows).to_csv(
        MANIFEST_DIR / "analysis_manifest_all.tsv", sep="\t", index=False
    )

    trait_columns = traits
    trait_aliases = [trait_alias[trait] for trait in trait_columns]
    exposure_aliases = [exposure_alias[exposure] for exposure in exposure_names]
    numeric_covariates = ["Age", "avMSE"] + [f"PC{index}" for index in range(1, 11)]

    keep = pd.read_csv(KEEP_FILE, sep=r"\s+", dtype={"IID": "string"})
    keep = keep.loc[:, ["IID"]].copy()
    keep["IID"] = normalize_iid(keep["IID"])

    phenotype = pd.read_csv(
        PHENO_FILE,
        sep=r"\s+",
        usecols=["IID"] + trait_columns,
        dtype={"IID": "string"},
        low_memory=False,
    )
    phenotype["IID"] = normalize_iid(phenotype["IID"])
    for column in trait_columns:
        phenotype[column] = pd.to_numeric(phenotype[column])

    covariate = pd.read_csv(
        COVAR_FILE,
        sep=r"\s+",
        usecols=["IID", "Sex", "gbatch"] + numeric_covariates,
        dtype={"IID": "string", "Sex": "string", "gbatch": "string"},
        low_memory=False,
    )
    covariate["IID"] = normalize_iid(covariate["IID"])
    covariate["Sex"] = pd.to_numeric(
        covariate["Sex"].str.strip().replace({"c0": 0, "c1": 1})
    )
    for column in numeric_covariates:
        covariate[column] = pd.to_numeric(covariate[column])

    exposure = pd.read_csv(
        EXPOSURE_FILE,
        usecols=["eid"] + exposure_names,
        dtype={"eid": "string"},
        low_memory=False,
    ).rename(columns={"eid": "IID"})
    exposure["IID"] = normalize_iid(exposure["IID"])
    for column in exposure_names:
        exposure[column] = pd.to_numeric(exposure[column])

    merged = (
        keep.merge(phenotype, on="IID")
        .merge(covariate, on="IID")
        .merge(exposure, on="IID")
        .dropna(subset=["Sex"] + numeric_covariates)
        .reset_index(drop=True)
    )
    batch = merged["gbatch"].astype("string").str.strip()
    batch_levels = sorted(batch.unique().tolist())
    batch_dummies = pd.get_dummies(batch, dtype="int8").reindex(
        columns=batch_levels[1:], fill_value=0
    )
    batch_columns = [f"B{index:03d}" for index in range(1, len(batch_levels))]
    batch_dummies.columns = batch_columns

    trait_data = merged.loc[:, trait_columns].copy()
    trait_data.columns = trait_aliases
    exposure_data = merged.loc[:, exposure_names].copy()
    exposure_data.columns = exposure_aliases
    covariate_data = merged.loc[
        :, ["IID", "Age", "Sex", "avMSE"] + [f"PC{index}" for index in range(1, 11)]
    ]
    master = pd.concat(
        [covariate_data.reset_index(drop=True), trait_data, exposure_data, batch_dummies],
        axis=1,
    )
    master.to_csv(INPUT_DIR / "discovery_gem_master.csv", index=False, na_rep="NA")

    covariate_names = (
        ["Age", "Sex", "avMSE"]
        + [f"PC{index}" for index in range(1, 11)]
        + batch_columns
    )
    (INPUT_DIR / "gem_covar_names.txt").write_text(
        "\n".join(covariate_names) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
