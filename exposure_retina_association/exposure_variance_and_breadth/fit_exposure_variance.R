# Incremental ordinary R-squared on identical complete cases; exposures retain original coding.
library(data.table)
ids <- data.table::fread("data/discovery_ids.tsv")
traits <- data.table::fread("data/representative_traits.tsv")
exposures <- data.table::fread("data/exposure_dictionary.csv")
exposures <- exposures[, .(exposure = association_exposure, expression_exposure)]
phenotype <- data.table::fread("data/discovery_retinal_phenotypes.tsv")
exposure <- data.table::fread("data/discovery_exposures.csv")
covariate <- data.table::fread("data/discovery_covariates.csv")
dat <- merge(ids, phenotype, by = "eid", all.x = TRUE, sort = FALSE)
dat <- merge(dat, exposure, by = "eid", all.x = TRUE, sort = FALSE)
dat <- merge(dat, covariate, by = "eid", all.x = TRUE, sort = FALSE)
pairs <- expand.grid(exposure = exposures$exposure, EyeFeature = traits$EyeFeature,
                     stringsAsFactors = FALSE)
result <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  trait <- pairs$EyeFeature[i]
  term <- pairs$exposure[i]
  column <- exposures$expression_exposure[match(term, exposures$exposure)]
  model_data <- data.frame(Y = dat[[trait]], X = dat[[column]], Age = dat$Age,
                          Sex = factor(dat$Sex), avMSE = dat$avMSE)
  model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]
  baseline <- stats::lm(Y ~ Age + Sex + avMSE, data = model_data, na.action = na.fail)
  full <- stats::lm(Y ~ Age + Sex + avMSE + X, data = model_data, na.action = na.fail)
  baseline_summary <- summary(baseline)
  full_summary <- summary(full)
  beta <- full_summary$coefficients["X", ]
  delta <- full_summary$r.squared - baseline_summary$r.squared
  result[[i]] <- data.table::data.table(exposure = term, EyeFeature = trait, n = stats::nobs(full),
    r2_baseline = baseline_summary$r.squared, r2_full = full_summary$r.squared,
    delta_r2 = delta, delta_r2_pct = 100 * delta,
    estimate = beta[1], std.error = beta[2], p.value = beta[4])
  if (i %% 100L == 0L || i == nrow(pairs)) {
    message(sprintf("Exposure variance: %d/%d (%.1f%%)", i, nrow(pairs), 100 * i / nrow(pairs)))
  }
}
result <- data.table::rbindlist(result)
result[, significant := p.value < 0.05 / (288 * 58)]
dir.create("results/exposure_variance", recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(result, "results/exposure_variance/pair_results.tsv", sep = "\t")
