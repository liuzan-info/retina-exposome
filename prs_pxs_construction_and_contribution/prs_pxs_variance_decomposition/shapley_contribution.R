# Two-score Shapley decomposition with Age and all other covariates held fixed.
library(dplyr)
library(tidyr)
metrics <- readr::read_tsv("results/prs_pxs_variance/model_metrics.tsv", show_col_types = FALSE)
dictionary <- readr::read_csv("data/retinal_dictionary.csv", show_col_types = FALSE)
wide <- metrics |>
  dplyr::select(EyeFeature, model, n, r2) |>
  tidyr::pivot_wider(names_from = model, values_from = c(n, r2))
stopifnot(all(wide$n_model_base == wide$n_model_base_pxs),
          all(wide$n_model_base == wide$n_model_base_prs),
          all(wide$n_model_base == wide$n_model_base_pxs_prs))
contribution <- wide |>
  dplyr::mutate(
    delta_r2_joint = r2_model_base_pxs_prs - r2_model_base,
    shapley_prs = 0.5 * ((r2_model_base_prs - r2_model_base) +
                         (r2_model_base_pxs_prs - r2_model_base_pxs)),
    shapley_pxs = 0.5 * ((r2_model_base_pxs - r2_model_base) +
                         (r2_model_base_pxs_prs - r2_model_base_prs)),
    prs_share = shapley_prs / delta_r2_joint,
    pxs_share = shapley_pxs / delta_r2_joint)
group_contribution <- contribution |>
  dplyr::left_join(dplyr::select(dictionary, EyeFeature, group), by = "EyeFeature") |>
  dplyr::group_by(group) |>
  dplyr::summarise(shapley_prs = mean(shapley_prs), shapley_pxs = mean(shapley_pxs), .groups = "drop") |>
  dplyr::mutate(prs_share = shapley_prs / (shapley_prs + shapley_pxs),
                pxs_share = shapley_pxs / (shapley_prs + shapley_pxs))
group_metrics <- metrics |>
  dplyr::left_join(dplyr::select(dictionary, EyeFeature, group), by = "EyeFeature") |>
  dplyr::group_by(group, model) |>
  dplyr::summarise(mean_r2 = mean(r2), .groups = "drop")
readr::write_tsv(contribution, "results/prs_pxs_variance/trait_shapley.tsv")
readr::write_tsv(group_contribution, "results/prs_pxs_variance/group_shapley.tsv")
readr::write_tsv(group_metrics, "results/prs_pxs_variance/group_model_metrics.tsv")
