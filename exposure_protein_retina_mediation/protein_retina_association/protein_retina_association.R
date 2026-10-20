# Discovery retinal phenotype ~ protein + Age + Sex + avMSE + eGFR.
# Inputs retain upstream phenotype normalization, imputed covariates and prepared protein names.
phenotypes <- readr::read_tsv("data/retinal_phenotypes.tsv", show_col_types = FALSE)
proteins <- readr::read_csv("data/discovery_proteins.csv", show_col_types = FALSE)
covariates <- readr::read_csv("data/discovery_protein_covariates.csv", show_col_types = FALSE)
ids <- read.table("data/discovery_eids.txt", header = FALSE, col.names = "IID")
traits <- readr::read_csv("data/retinal_dictionary.csv", show_col_types = FALSE)$EyeFeature
protein_names <- setdiff(names(proteins), "eid")
covars <- c("Age", "Sex", "avMSE", "eGFR")
data <- phenotypes |>
  dplyr::filter(IID %in% ids$IID) |>
  dplyr::select(IID, dplyr::all_of(traits)) |>
  dplyr::inner_join(dplyr::select(covariates, eid, dplyr::all_of(covars)), by = c("IID" = "eid")) |>
  dplyr::inner_join(proteins, by = c("IID" = "eid")) |>
  dplyr::mutate(Sex = factor(Sex))
results <- list()
for (i in seq_along(traits)) {
  trait <- traits[i]
  current <- data[!is.na(data[[trait]]), , drop = FALSE]
  valid_covars <- covars[vapply(current[covars], function(x) length(unique(stats::na.omit(x))) >= 2L, logical(1))]
  for (protein in protein_names) {
    values <- stats::na.omit(current[[protein]])
    if (length(unique(values)) < 2L) next
    if (is.numeric(values) && stats::sd(values) == 0) next
    formula <- stats::as.formula(paste0("`", trait, "` ~ ",
      paste0("`", c(valid_covars, protein), "`", collapse = " + ")))
    fit <- stats::lm(formula, data = current)
    coefficient <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95) |>
      dplyr::mutate(term = gsub("`", "", term, fixed = TRUE)) |>
      dplyr::filter(term == protein) |>
      dplyr::mutate(protein = protein, EyeFeature = trait, nobs = stats::nobs(fit))
    results[[length(results) + 1L]] <- coefficient
  }
  message(sprintf("Protein-retinal models: %d/%d (%.1f%%)", i, length(traits), 100 * i / length(traits)))
}
results <- dplyr::bind_rows(results) |>
  dplyr::mutate(p.value.fdr = stats::p.adjust(p.value, "BH"),
                p.value.bonferroni = stats::p.adjust(p.value, "bonferroni"))
dir.create("results/protein_retina", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(results, "results/protein_retina/protein_retina_association.csv")
