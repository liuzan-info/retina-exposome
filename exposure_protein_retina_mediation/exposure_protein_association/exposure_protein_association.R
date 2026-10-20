# Protein ~ exposure + Age + Sex + avMSE + eGFR in processed discovery data.
proteins <- readr::read_csv("data/discovery_proteins.csv", show_col_types = FALSE)
exposures <- readr::read_csv("data/discovery_exposures.csv", show_col_types = FALSE)
covariates <- readr::read_csv("data/discovery_protein_covariates.csv", show_col_types = FALSE)
dictionary <- readr::read_csv("data/exposure_dictionary.csv", show_col_types = FALSE)
covars <- c("Age", "Sex", "avMSE", "eGFR")
protein_names <- setdiff(names(proteins), "eid")
data <- proteins |>
  dplyr::left_join(dplyr::select(covariates, eid, dplyr::all_of(covars)), by = "eid") |>
  dplyr::left_join(dplyr::select(exposures, eid, dplyr::all_of(dictionary$expression_exposure)), by = "eid") |>
  dplyr::mutate(Sex = factor(Sex))
results <- list()
for (i in seq_along(protein_names)) {
  protein <- protein_names[i]
  current <- data[!is.na(data[[protein]]), , drop = FALSE]
  valid_covars <- covars[vapply(current[covars], function(x) length(unique(stats::na.omit(x))) >= 2L, logical(1))]
  for (j in seq_len(nrow(dictionary))) {
    exposure <- dictionary$expression_exposure[j]
    values <- stats::na.omit(current[[exposure]])
    if (length(unique(values)) < 2L) next
    if (is.numeric(values) && stats::sd(values) == 0) next
    formula <- stats::as.formula(paste0("`", protein, "` ~ ",
      paste0("`", c(valid_covars, exposure), "`", collapse = " + ")))
    fit <- stats::lm(formula, data = current)
    coefficient <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95) |>
      dplyr::mutate(term = gsub("`", "", term, fixed = TRUE)) |>
      dplyr::filter(term == exposure) |>
      dplyr::mutate(term = dictionary$association_exposure[j],
        exposure = dictionary$association_exposure[j], protein = protein, nobs = stats::nobs(fit))
    results[[length(results) + 1L]] <- coefficient
  }
  if (i %% 20L == 0L || i == length(protein_names)) {
    message(sprintf("Protein-exposure models: %d/%d (%.1f%%)", i, length(protein_names), 100 * i / length(protein_names)))
  }
}
results <- dplyr::bind_rows(results) |>
  dplyr::mutate(p.value.fdr = stats::p.adjust(p.value, "BH"),
                p.value.bonferroni = stats::p.adjust(p.value, "bonferroni"))
dir.create("results/exposure_protein", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(results, "results/exposure_protein/protein_exposure_association.csv")
