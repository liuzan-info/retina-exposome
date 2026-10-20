#!/usr/bin/env python3
# No-significant-clump chromosomes produce no .clumps file, as recorded by PLINK.
from pathlib import Path
import pandas as pd
root = Path("results/vqtl")
rows = []
for trait in (root / "pheno/phenotype_list.txt").read_text().splitlines():
    for chrom in range(1, 23):
        prefix = root / "clumps" / trait / f"chr{chrom}"
        if "No significant --clump results.  Skipping." in prefix.with_suffix(".log").read_text():
            continue
        frame = pd.read_csv(prefix.with_suffix(".clumps"), sep="\t")
        frame.insert(0, "safe_phenotype", trait)
        rows.append(frame)
pd.concat(rows, ignore_index=True).to_csv(root / "all_58phenotypes.clumps", sep="\t", index=False)
