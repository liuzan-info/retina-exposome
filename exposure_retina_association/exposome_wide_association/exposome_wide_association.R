# Discovery exposome-wide association, restricted to the final 288 exposures and 58 phenotypes.
# Prepared phenotypes/exposures retain the original dataset-specific normalization and coding.
phenotypes <- readr::read_tsv("data/retinal_phenotypes.tsv", show_col_types = FALSE)
exposures <- readr::read_csv("data/discovery_exposures.csv", show_col_types = FALSE)
covariates <- readr::read_csv("data/discovery_covariates.csv", show_col_types = FALSE)
ids <- read.table("data/discovery_eids.txt", header = FALSE, col.names = "IID")
dictionary <- readr::read_csv("data/exposure_dictionary.csv", show_col_types = FALSE)
traits <- readr::read_csv("data/retinal_dictionary.csv", show_col_types = FALSE)$EyeFeature
covars <- c("Age", "Sex", "avMSE")
data <- phenotypes |>
  dplyr::filter(IID %in% ids$IID) |>
  dplyr::select(IID, dplyr::all_of(traits)) |>
  dplyr::left_join(dplyr::select(covariates, eid, dplyr::all_of(covars)), by = c("IID" = "eid")) |>
  dplyr::left_join(dplyr::select(exposures, eid, dplyr::all_of(dictionary$expression_exposure)), by = c("IID" = "eid")) |>
  dplyr::mutate(Sex = factor(Sex))
results <- list()
for (i in seq_along(traits)) {
  trait <- traits[i]
  current <- data[!is.na(data[[trait]]), , drop = FALSE]
  valid_covars <- covars[vapply(current[covars], function(x) length(unique(stats::na.omit(x))) >= 2L, logical(1))]
  for (j in seq_len(nrow(dictionary))) {
    exposure <- dictionary$expression_exposure[j]
    values <- stats::na.omit(current[[exposure]])
    if (length(unique(values)) < 2L) next
    if (is.numeric(values) && stats::sd(values) == 0) next
    formula <- stats::as.formula(paste0("`", trait, "` ~ ",
      paste0("`", c(valid_covars, exposure), "`", collapse = " + ")))
    fit <- stats::lm(formula, data = current)
    coefficient <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95) |>
      dplyr::mutate(term = gsub("`", "", term, fixed = TRUE)) |>
      dplyr::filter(term == exposure) |>
      dplyr::mutate(term = dictionary$association_exposure[j],
        oct_features = trait, nobs = stats::nobs(fit))
    results[[length(results) + 1L]] <- coefficient
  }
  message(sprintf("Exposome-wide models: %d/%d (%.1f%%)", i, length(traits), 100 * i / length(traits)))
}
results <- dplyr::bind_rows(results) |>
  dplyr::mutate(p.value.bonferroni = pmin(p.value * (288 * 58), 1),
                p.value.fdr = stats::p.adjust(p.value, "BH"),
                significant = p.value < 0.05 / (288 * 58))
dir.create("results/exposome_wide_association", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(results, "results/exposome_wide_association/discovery_exposure_associations.csv")
