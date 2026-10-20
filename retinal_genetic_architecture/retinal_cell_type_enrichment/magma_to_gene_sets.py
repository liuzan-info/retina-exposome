#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def load_gene_map(path):
    table = pd.read_csv(path, sep=r"\s+", header=None)
    return table[[0, table.columns[-1]]].rename(columns={0: "GENE", table.columns[-1]: "SYMBOL"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--gene-loc", required=True)
    parser.add_argument("--magma-gene-dir", required=True)
    parser.add_argument("--zscore-out", required=True)
    args = parser.parse_args()

    traits = pd.read_csv(args.manifest, sep="\t")["trait"].tolist()
    gene_map = load_gene_map(args.gene_loc)
    zscores = []

    for trait in traits:
        gene_file = Path(args.magma_gene_dir) / f"{trait}.genes.out"
        result = pd.read_csv(gene_file, sep=r"\s+")[["GENE", "ZSTAT"]]
        result = result.merge(gene_map, on="GENE", how="inner")
        result = result.drop_duplicates("SYMBOL").set_index("SYMBOL")[["ZSTAT"]]
        result.columns = [trait]
        zscores.append(result)

    wide = pd.concat(zscores, axis=1, join="inner").reset_index()
    wide = wide.rename(columns={"SYMBOL": "GENE"})
    wide.to_csv(args.zscore_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
