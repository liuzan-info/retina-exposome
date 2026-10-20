# Discovery-weighted total and category PXS in the aligned replication population.
# Exposure values and discovery_estimate retain the original statistical coding.
library(data.table)
exposure <- data.table::fread("data/scores/01_rep1_exposure.tsv")
mapping <- data.table::fread("data/selected_exposures.csv")
out <- "results/pxs"
dir.create(file.path(out, "group_z"), recursive = TRUE, showWarnings = FALSE)
zscore_strict <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  x_sd <- sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / x_sd
}


compute_weighted_score <- function(exposure_dt, exposure_cols, betas) {
  exposure_matrix <- as.matrix(exposure_dt[, ..exposure_cols])
  contribution_matrix <- sweep(exposure_matrix, 2L, betas, `*`)
  non_missing_n <- rowSums(!is.na(exposure_matrix))
  weighted_score <- rowSums(contribution_matrix, na.rm = TRUE) / non_missing_n
  weighted_score[non_missing_n == 0] <- NA_real_
  as.numeric(weighted_score)
}

traits <- unique(mapping$EyeFeature)
total_raw <- data.table::data.table(IID = exposure$IID)
total_z <- data.table::data.table(IID = exposure$IID)
for (i in seq_along(traits)) {
  trait <- traits[i]
  selected <- mapping[EyeFeature == trait]
  values <- compute_weighted_score(exposure, selected$expression_exposure, selected$discovery_estimate)
  total_raw[[trait]] <- values
  total_z[[trait]] <- zscore_strict(values)
  message(sprintf("PXS: %d/%d (%.1f%%)", i, length(traits), 100 * i / length(traits)))
}
data.table::fwrite(total_raw, file.path(out, "03_total_pxs_raw.tsv"), sep = "\t")
data.table::fwrite(total_z, file.path(out, "03_total_pxs_z.tsv"), sep = "\t")
categories <- unique(mapping$Category)
manifest <- vector("list", length(categories))
for (i in seq_along(categories)) {
  category <- categories[i]
  score <- data.table::data.table(IID = exposure$IID)
  for (trait in traits) {
    selected <- mapping[EyeFeature == trait & Category == category]
    score[[trait]] <- if (nrow(selected) == 0L) rep(NA_real_, nrow(exposure)) else {
      zscore_strict(compute_weighted_score(exposure, selected$expression_exposure, selected$discovery_estimate))
    }
  }
  key <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(category)))
  path <- file.path(out, "group_z", paste0(key, ".tsv"))
  data.table::fwrite(score, path, sep = "\t")
  manifest[[i]] <- data.table::data.table(group = category, path = path)
}
data.table::fwrite(data.table::rbindlist(manifest), file.path(out, "03_group_pxs_z.tsv"), sep = "\t")
