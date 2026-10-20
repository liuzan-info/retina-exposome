#!/usr/bin/env python3

import argparse
import gc
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from anndata._core.sparse_dataset import SparseDataset


def major_cell_type(cell_type):
    label = cell_type.lower()
    if "pigment epithelial" in label:
        return "RPE"
    if any(value in label for value in ["mueller", "muller", "müller"]):
        return "Muller_glia"
    if "microglia" in label:
        return "microglia"
    if "astrocyte" in label:
        return "astrocyte"
    if "ganglion" in label:
        return "RGC"
    if "offx cell" in label:
        return "bipolar"
    if "bipolar" in label:
        return "bipolar"
    if "amacrine" in label:
        return "amacrine"
    if "horizontal" in label:
        return "horizontal"
    if "cone" in label:
        return "cone"
    if "rod" in label:
        return "rod"
    raise ValueError(f"Unclassified HRCA cell type: {cell_type}")


def region_label(tissue):
    if tissue == "fovea centralis":
        return "fovea"
    if tissue in {"macula lutea", "macula lutea proper"}:
        return "macula"
    if tissue == "peripheral region of retina":
        return "peripheral"
    raise ValueError(f"Unclassified HRCA tissue: {tissue}")


def is_adult(stage):
    label = stage.lower()
    match = re.search(r"(\d+)[ -]year-old", label)
    if match is not None:
        return int(match.group(1)) >= 18
    return any(value in label for value in ["adult", "decade", "middle aged"])


def write_raw_subset(input_path, output_path, source_positions, raw_gene_keep):
    n_obs = len(source_positions)
    n_var = int(raw_gene_keep.sum())
    n_genes = np.empty(n_obs, dtype=np.int32)
    n_counts = np.empty(n_obs, dtype=np.float64)
    indptr = np.empty(n_obs + 1, dtype=np.int32)
    indptr[0] = 0

    with h5py.File(input_path, "r") as source_file, h5py.File(
        output_path, "r+"
    ) as output_file:
        raw_x = SparseDataset(source_file["raw/X"])
        del output_file["X"]
        output_x = output_file.create_group("X")
        output_x.attrs["encoding-type"] = "csr_matrix"
        output_x.attrs["encoding-version"] = "0.1.0"
        output_x.attrs["shape"] = np.array([n_obs, n_var], dtype=np.int64)
        output_data = output_x.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            chunks=(1_000_000,),
            dtype=np.float32,
            compression="gzip",
        )
        output_indices = output_x.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(None,),
            chunks=(1_000_000,),
            dtype=np.int32,
            compression="gzip",
        )

        value_offset = 0
        row_offset = 0
        for start in range(0, n_obs, 5000):
            chunk = raw_x[source_positions[start : start + 5000], :][
                :, raw_gene_keep
            ].tocsr()
            chunk.sort_indices()
            next_value_offset = value_offset + chunk.nnz
            output_data.resize((next_value_offset,))
            output_indices.resize((next_value_offset,))
            output_data[value_offset:next_value_offset] = chunk.data
            output_indices[value_offset:next_value_offset] = chunk.indices

            next_row_offset = row_offset + chunk.shape[0]
            indptr[row_offset + 1 : next_row_offset + 1] = (
                value_offset + chunk.indptr[1:]
            )
            n_genes[row_offset:next_row_offset] = np.diff(chunk.indptr)
            n_counts[row_offset:next_row_offset] = np.asarray(
                chunk.sum(axis=1)
            ).ravel()
            value_offset = next_value_offset
            row_offset = next_row_offset
            print(
                f"rows {row_offset:,}/{n_obs:,}: {value_offset:,} nonzero values",
                flush=True,
            )

        output_x.create_dataset("indptr", data=indptr, compression="gzip")

    return n_genes, n_counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-per-stratum", type=int, default=200)
    parser.add_argument("--seeds", default="1,2,4,5")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    source = ad.read_h5ad(args.input, backed="r")
    obs = source.obs.copy()
    obs["_source_index"] = np.arange(source.n_obs)
    obs = obs.loc[obs["disease"].eq("normal")]
    obs = obs.loc[obs["development_stage"].astype(str).map(is_adult)]
    obs["major_cell_type"] = obs["cell_type"].astype(str).map(major_cell_type)
    obs["cell_subtype"] = obs["cell_type"].astype(str)
    obs["region"] = obs["tissue"].astype(str).map(region_label)
    obs["donor"] = obs["donor_id"].astype(str)
    obs["major_cell_subtype"] = obs["major_cell_type"] + "|" + obs["cell_subtype"]

    strata = list(
        obs.groupby(
            ["donor", "region", "major_cell_type"], observed=True, sort=True
        ).indices.values()
    )
    raw_gene_names = source.raw.var["feature_name"].astype(str)
    raw_gene_keep = (~raw_gene_names.duplicated()).to_numpy()
    raw_var = source.raw.var.loc[raw_gene_keep].copy()
    raw_var.index = raw_gene_names.loc[raw_gene_keep].to_numpy()
    raw_var.index.name = None

    selections = {}
    for seed in [int(value) for value in args.seeds.split(",")]:
        rng = np.random.default_rng(seed)
        selected_positions = np.sort(
            np.concatenate(
                [
                    rng.choice(
                        positions,
                        size=min(len(positions), args.max_per_stratum),
                        replace=False,
                    )
                    for positions in strata
                ]
            )
        )
        selected_obs = obs.iloc[selected_positions].copy()
        source_positions = selected_obs.pop("_source_index").to_numpy()
        selected_obs = selected_obs[
            [
                "major_cell_type",
                "cell_subtype",
                "region",
                "donor",
                "major_cell_subtype",
            ]
        ]
        partial_path = out_dir / f".hrca_seed{seed}.h5ad.partial"
        template = ad.AnnData(
            X=sparse.csr_matrix(
                (len(source_positions), int(raw_gene_keep.sum())), dtype=np.float32
            ),
            obs=selected_obs,
            var=raw_var.copy(),
        )
        template.write_h5ad(partial_path, compression="gzip")
        selections[seed] = (source_positions, partial_path)
        print(f"seed {seed}: selected {len(source_positions):,} cells", flush=True)

    source.file.close()
    del source, obs, strata, template, selected_obs
    gc.collect()

    for seed, (source_positions, partial_path) in selections.items():
        n_genes, n_counts = write_raw_subset(
            args.input, partial_path, source_positions, raw_gene_keep
        )
        prepared = ad.read_h5ad(partial_path, backed="r")
        obs_names = prepared.obs_names.copy()
        prepared.file.close()
        covariates = pd.DataFrame(
            {
                "const": 1.0,
                "n_genes": n_genes,
                "log1p_n_counts": np.log1p(n_counts),
            },
            index=obs_names,
        )
        covariates.to_csv(out_dir / f"hrca_seed{seed}.cov", sep="\t", index_label="index")
        partial_path.replace(out_dir / f"hrca_seed{seed}.h5ad")
        print(f"seed {seed}: files written", flush=True)


if __name__ == "__main__":
    main()
