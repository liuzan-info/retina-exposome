# Apply separate BH corrections across all 4,528 candidate paths before filtering.
raw <- readr::read_csv("results/protein_mediation/mediation_results_boot2000.csv", show_col_types = FALSE)
results <- raw |>
  dplyr::mutate(acme_fdr_bh = stats::p.adjust(acme_p, method = "BH"),
                total_fdr_bh = stats::p.adjust(total_p, method = "BH"))
significant <- results |> dplyr::filter(acme_fdr_bh < 0.05, total_fdr_bh < 0.05)
protein_counts <- significant |>
  dplyr::count(protein, name = "n_paths") |>
  dplyr::arrange(dplyr::desc(n_paths), protein)
readr::write_csv(results, "results/protein_mediation/mediation_results_adjusted.csv")
readr::write_csv(significant, "results/protein_mediation/significant_mediation_paths.csv")
readr::write_csv(protein_counts, "results/protein_mediation/protein_path_counts.csv")
