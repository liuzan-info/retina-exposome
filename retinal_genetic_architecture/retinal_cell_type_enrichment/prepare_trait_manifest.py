#!/usr/bin/env python3

import argparse
import pandas as pd


CATEGORY_BY_GROUP = {
    "Macular": "neural_retinal_thickness",
    "RNFL": "neural_retinal_thickness",
    "GC-IPL": "neural_retinal_thickness",
    "INL": "neural_retinal_thickness",
    "INL-ELM": "neural_retinal_thickness",
    "ELM-ISOS": "neural_retinal_thickness",
    "ISOS-RPE": "neural_retinal_thickness",
    "INL-RPE": "neural_retinal_thickness",
    "RPE": "rpe_thickness",
    "Vessel caliber": "vascular",
    "Vascular tortuosity": "vascular",
    "Vascular fractal dimension": "vascular",
    "Papillomacular arcade": "vascular",
    "Cup-to-disc ratio": "optic_disc",
    "Optic disc morphology": "optic_disc",
}

SUITABILITY = {
    "neural_retinal_thickness": "high",
    "rpe_thickness": "expression_match_power_limited",
    "vascular": "insufficient_without_vascular_atlas",
    "optic_disc": "partial_retinal_component_only",
}

EXPECTED_CELLS_BY_GROUP = {
    "Macular": "rod;cone;Muller_glia;RGC;RPE",
    "RNFL": "RGC;Muller_glia;astrocyte",
    "GC-IPL": "RGC;bipolar;amacrine;Muller_glia",
    "INL": "bipolar;amacrine;horizontal;Muller_glia",
    "INL-ELM": "bipolar;amacrine;horizontal;Muller_glia",
    "ELM-ISOS": "rod;cone;Muller_glia",
    "ISOS-RPE": "rod;cone;RPE;Muller_glia",
    "INL-RPE": "rod;cone;bipolar;Muller_glia;RPE",
    "RPE": "RPE",
    "Vessel caliber": "endothelial;pericyte;vascular_smooth_muscle",
    "Vascular tortuosity": "endothelial;pericyte;vascular_smooth_muscle",
    "Vascular fractal dimension": "endothelial;pericyte;vascular_smooth_muscle",
    "Papillomacular arcade": "endothelial;pericyte;vascular_smooth_muscle",
    "Cup-to-disc ratio": "RGC;astrocyte;Muller_glia",
    "Optic disc morphology": "RGC;astrocyte;Muller_glia",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-file", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    traits = pd.read_csv(args.trait_file, encoding="utf-8-sig")
    traits["trait"] = traits["EyeFeature"]
    traits["category"] = traits["group"].map(CATEGORY_BY_GROUP)
    traits["h2_z"] = traits["ldsc_h2_obs_over_se_with_merge_alleles"]
    traits["analysis_tier"] = traits["h2_z"].ge(5).map(
        {True: "validation", False: "exploratory"}
    )
    traits["hrca_suitability"] = traits["category"].map(SUITABILITY)
    traits["expected_cell_types"] = traits["group"].map(EXPECTED_CELLS_BY_GROUP)
    traits["hrca_has_expected_vascular_cells"] = pd.Series(
        pd.NA, index=traits.index, dtype="boolean"
    )
    traits.loc[
        traits["category"].eq("vascular"), "hrca_has_expected_vascular_cells"
    ] = False

    columns = [
        "trait",
        "feature_name",
        "plot_label_med",
        "group",
        "category",
        "h2_z",
        "analysis_tier",
        "hrca_suitability",
        "expected_cell_types",
        "hrca_has_expected_vascular_cells",
    ]
    traits[columns].to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
