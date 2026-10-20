# Four nested score models on the same complete cases for each phenotype.
# Age remains a fixed covariate; no three-variable Shapley decomposition is performed.
feature <- read.delim("data/scores/01_rep1_oct_feature.tsv", check.names = FALSE)
covariate <- read.delim("data/scores/01_rep1_covariate.tsv", check.names = FALSE)
prs <- read.delim("data/scores/01_rep1_prs_avg.tsv", check.names = FALSE)
pxs <- read.delim("data/scores/03_total_pxs_z.tsv", check.names = FALSE)
stopifnot(identical(feature$IID, covariate$IID), identical(feature$IID, prs$IID),
          identical(feature$IID, pxs$IID))
covariate$Age <- as.numeric(scale(covariate$Age))
covariate$avMSE <- as.numeric(scale(covariate$avMSE))
covariate$Sex <- factor(covariate$Sex)
covariate$gbatch <- factor(covariate$gbatch)
base_terms <- c("Age", "Sex", "avMSE", "gbatch", paste0("PC", 1:10))
for (pc in paste0("PC", 1:10)) covariate[[pc]] <- as.numeric(covariate[[pc]])
traits <- intersect(setdiff(names(feature), "IID"), setdiff(names(pxs), "IID"))
model_terms <- list(model_base = character(), model_base_pxs = "PXS",
                    model_base_prs = "PRS", model_base_pxs_prs = c("PRS", "PXS"))
metrics <- comparisons <- vector("list", length(traits))
for (i in seq_along(traits)) {
  trait <- traits[i]
  dat <- data.frame(Y = feature[[trait]], covariate[, base_terms],
                    PRS = prs[[paste0(trait, "_AVG")]], PXS = pxs[[trait]])
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  fits <- lapply(model_terms, function(terms) {
    stats::lm(stats::reformulate(c(base_terms, terms), response = "Y"),
              data = dat, na.action = na.fail)
  })
  r2 <- vapply(fits, function(fit) summary(fit)$r.squared, numeric(1))
  metrics[[i]] <- data.frame(EyeFeature = trait, model = names(fits),
                            n = nrow(dat), r2 = unname(r2))
  comparisons[[i]] <- data.frame(
    EyeFeature = trait, n = nrow(dat),
    delta_r2_prs_given_pxs = r2["model_base_pxs_prs"] - r2["model_base_pxs"],
    delta_r2_pxs_given_prs = r2["model_base_pxs_prs"] - r2["model_base_prs"],
    delta_r2_joint_vs_base = r2["model_base_pxs_prs"] - r2["model_base"],
    p_prs_given_pxs = stats::anova(fits$model_base_pxs, fits$model_base_pxs_prs)$`Pr(>F)`[2],
    p_pxs_given_prs = stats::anova(fits$model_base_prs, fits$model_base_pxs_prs)$`Pr(>F)`[2],
    p_joint_vs_base = stats::anova(fits$model_base, fits$model_base_pxs_prs)$`Pr(>F)`[2],
    row.names = NULL)
  message(sprintf("Score models: %d/%d (%.1f%%)", i, length(traits), 100 * i / length(traits)))
}
metrics <- do.call(rbind, metrics)
comparisons <- do.call(rbind, comparisons)
comparisons$fdr_prs_given_pxs <- stats::p.adjust(comparisons$p_prs_given_pxs, "BH")
comparisons$fdr_pxs_given_prs <- stats::p.adjust(comparisons$p_pxs_given_prs, "BH")
comparisons$fdr_joint_vs_base <- stats::p.adjust(comparisons$p_joint_vs_base, "BH")
dir.create("results/prs_pxs_variance", recursive = TRUE, showWarnings = FALSE)
write.table(metrics, "results/prs_pxs_variance/model_metrics.tsv", sep = "\t", row.names = FALSE)
write.table(comparisons, "results/prs_pxs_variance/model_comparisons.tsv", sep = "\t", row.names = FALSE)
