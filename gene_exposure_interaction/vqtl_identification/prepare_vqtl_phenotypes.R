# Raw phenotype residualization without Sex; remove |z|>5, then rescale retained residuals within Sex.
raw <- read.table("data/raw_retinal_phenotypes.tsv", header = TRUE, check.names = FALSE)
covariate <- read.table("data/gwas_covariates.tsv", header = TRUE, check.names = FALSE)
discovery <- read.table("data/discovery_keep.tsv", header = TRUE, check.names = FALSE)
traits <- read.csv("data/retinal_dictionary.csv", check.names = FALSE)$EyeFeature
covariate$IID <- as.character(covariate$IID)
raw$IID <- as.character(raw$IID)
discovery$IID <- as.character(discovery$IID)
base <- covariate[match(discovery$IID, covariate$IID), , drop = FALSE]
raw_index <- match(base$IID, raw$IID)
base$Sex <- as.character(base$Sex)
stopifnot(identical(sort(unique(stats::na.omit(base$Sex))), c("c0", "c1")))
model_formula <- stats::as.formula(paste("raw_value ~ Age + avMSE + factor(gbatch) +",
                                       paste0("PC", 1:10, collapse = " + ")))
dir.create("results/vqtl/pheno", recursive = TRUE, showWarnings = FALSE)
safe_traits <- gsub("[^A-Za-z0-9._-]+", "_", traits)
safe_traits <- gsub("^[._-]+|[._-]+$", "", safe_traits)
for (i in seq_along(traits)) {
  dat <- base
  dat$raw_value <- raw[[traits[i]]][raw_index]
  columns <- c("raw_value", "Sex", "Age", "avMSE", "gbatch", paste0("PC", 1:10))
  dat <- dat[stats::complete.cases(dat[, columns, drop = FALSE]), , drop = FALSE]
  fit <- stats::lm(model_formula, data = dat)
  dat$residual <- stats::residuals(fit)
  dat$residual_z <- (dat$residual - mean(dat$residual)) / stats::sd(dat$residual)
  dat <- dat[abs(dat$residual_z) <= 5, , drop = FALSE]
  dat$PHENO <- NA_real_
  for (sex in c("c0", "c1")) {
    index <- dat$Sex == sex
    residual <- dat$residual[index]
    dat$PHENO[index] <- (residual - mean(residual)) / stats::sd(residual)
  }
  write.table(dat[, c("FID", "IID", "PHENO")],
    file.path("results/vqtl/pheno", paste0(safe_traits[i], ".pheno")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(sprintf("vQTL phenotypes: %d/%d (%.1f%%)", i, length(traits), 100 * i / length(traits)))
}
writeLines(safe_traits, "results/vqtl/pheno/phenotype_list.txt")
