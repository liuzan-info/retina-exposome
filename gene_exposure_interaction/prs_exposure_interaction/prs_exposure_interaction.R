suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

exposure_dictionary_path <- "data/exposure_dictionary.csv"
oct_dictionary_path <- "data/retinal_dictionary.csv"
discovery_association_path <- "data/discovery_exposure_associations.csv"
rep1_exposure_path <- "data/rep1_exposures.csv"
oct_feature_path <- "data/retinal_phenotypes.tsv"
covariate_path <- "data/gwas_covariates.tsv"
prs_path <- "data/rep1_prs.tsv"
result_dir <- "results/prs_exposure_interaction"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- 40L
primary_definition_value <- "discovery_main_effect_bonferroni"

interaction_with_main_effect_p_path <- file.path(result_dir, "01_interaction_with_prs_exposure_p.tsv")
main_effect_filtered_path <- file.path(result_dir, "02_prs_exposure_p_lt_0.05.tsv")
main_effect_filter_summary_path <- file.path(result_dir, "03_main_effect_filter_summary.tsv")
top_filtered_interactions_path <- file.path(result_dir, "04_top_filtered_interactions.tsv")

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.table(x), path, sep = "\t", na = "NA")
}

ordered_unique <- function(x) {
  x[!duplicated(x)]
}

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

align_to_ids <- function(dt, ids, id_col = "IID") {
  matched_index <- match(ids, dt[[id_col]])
  aligned_dt <- dt[matched_index]
  aligned_dt[[id_col]] <- ids
  aligned_dt
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  names(sort(table(x), decreasing = TRUE))[1]
}

safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
}

collapse_level_counts <- function(x, max_levels = 20L) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  tab <- table(x)
  if (length(tab) > max_levels) {
    return(paste0("n_unique_gt_", max_levels))
  }
  paste(paste0(names(tab), ":", as.integer(tab)), collapse = ";")
}

count_min_group_n <- function(x, max_levels = 20L) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_integer_)
  }
  tab <- table(x)
  if (length(tab) > max_levels) {
    return(NA_integer_)
  }
  as.integer(min(tab))
}

process_exposure_vector <- function(
  raw_vec,
  exposure_name,
  exposure_type,
  semantic_numeric_opposite
) {
  raw_num <- as.numeric(raw_vec)

  if (exposure_name == "e1210") {
    model_vec <- ifelse(raw_num == 1, 1, ifelse(raw_num == 2, 0, NA_real_))
    return(list(
      model_vec = model_vec,
      exposure_processed = "binary_recoded",
      exposure_original_coding = "1 = Yes; 2 = No",
      exposure_model_coding = "1 -> 1; 2 -> 0",
      semantic_transformation = "e1210_yes_1_no_0_recode",
      reference_level = "2 = No",
      exposed_level = "1 = Yes"
    ))
  }

  if (exposure_type == "CatOrd") {
    if (isTRUE(semantic_numeric_opposite)) {
      model_vec <- min(raw_num, na.rm = TRUE) + max(raw_num, na.rm = TRUE) - raw_num
      return(list(
        model_vec = zscore_strict(model_vec),
        exposure_processed = "z_scaled_ordered_numeric_semantic_mirrored",
        exposure_original_coding = "ordered numeric coding from exposure matrix",
        exposure_model_coding = "min(x) + max(x) - x, then z-scaled ordered numeric",
        semantic_transformation = "min(x) + max(x) - x",
        reference_level = NA_character_,
        exposed_level = NA_character_
      ))
    }

    return(list(
      model_vec = zscore_strict(raw_num),
      exposure_processed = "z_scaled_ordered_numeric",
      exposure_original_coding = "ordered numeric coding from exposure matrix",
      exposure_model_coding = "z-scaled ordered numeric",
      semantic_transformation = "none",
      reference_level = NA_character_,
      exposed_level = NA_character_
    ))
  }

  if (exposure_type %in% c("Binary", "CatUnord")) {
    return(list(
      model_vec = raw_num,
      exposure_processed = "raw_binary_or_dummy",
      exposure_original_coding = "0/1 coding from exposure matrix",
      exposure_model_coding = "0/1",
      semantic_transformation = "none",
      reference_level = "0",
      exposed_level = "1"
    ))
  }

  list(
    model_vec = raw_num,
    exposure_processed = "preprocessed_continuous",
    exposure_original_coding = "continuous exposure matrix value",
    exposure_model_coding = "preprocessed continuous value",
    semantic_transformation = "none",
    reference_level = NA_character_,
    exposed_level = NA_character_
  )
}

exposure_dictionary_dt <- fread(exposure_dictionary_path, check.names = FALSE)
oct_dictionary_dt <- fread(oct_dictionary_path, check.names = FALSE)
setnames(
  oct_dictionary_dt,
  sub(paste0("^", intToUtf8(0xfeff)), "", names(oct_dictionary_dt))
)

eye_features <- ordered_unique(oct_dictionary_dt$EyeFeature)
step3_exposures <- exposure_dictionary_dt[, .(
  expression_exposure,
  association_exposure,
  exposure_modifiability,
  Category,
  Type,
  semantic_numeric_opposite,
  FieldID,
  Field,
  exposure_value_meaning_simple,
  exposure_value_meaning
)]
stopifnot(
  nrow(step3_exposures) == 288L,
  uniqueN(step3_exposures$association_exposure) == 288L,
  length(eye_features) == 58L
)

prs_raw_dt <- fread(prs_path, check.names = FALSE)
prs_raw_dt[, IID := as.character(IID)]
analysis_ids <- ordered_unique(prs_raw_dt$IID[nzchar(prs_raw_dt$IID)])

oct_feature_raw_dt <- fread(oct_feature_path, check.names = FALSE)
oct_feature_raw_dt[, IID := as.character(IID)]

covariate_raw_dt <- fread(covariate_path, check.names = FALSE)
covariate_raw_dt[, IID := as.character(IID)]

exposure_raw_dt <- fread(rep1_exposure_path, check.names = FALSE)
exposure_raw_dt[, eid := as.character(eid)]
setnames(exposure_raw_dt, "eid", "IID")

prs_avg_cols <- paste0(eye_features, "_AVG")
missing_prs_cols <- setdiff(prs_avg_cols, names(prs_raw_dt))
if (length(missing_prs_cols) > 0) {
  stop("Missing PRS AVG columns: ", paste(missing_prs_cols, collapse = ", "))
}

missing_eye_feature_cols <- setdiff(eye_features, names(oct_feature_raw_dt))
if (length(missing_eye_feature_cols) > 0) {
  stop("Missing OCT feature columns: ", paste(missing_eye_feature_cols, collapse = ", "))
}

covariate_cols <- c("Age", "Sex", "avMSE", "gbatch", paste0("PC", 1:10))
missing_covariate_cols <- setdiff(c("IID", covariate_cols), names(covariate_raw_dt))
if (length(missing_covariate_cols) > 0) {
  stop("Missing covariate columns: ", paste(missing_covariate_cols, collapse = ", "))
}

present_exposure_cols <- intersect(step3_exposures$expression_exposure, names(exposure_raw_dt))
missing_exposure_cols <- setdiff(step3_exposures$expression_exposure, names(exposure_raw_dt))

feature_dt <- align_to_ids(
  oct_feature_raw_dt[, c("IID", eye_features), with = FALSE],
  analysis_ids,
  "IID"
)
covariate_dt <- align_to_ids(
  covariate_raw_dt[, c("IID", covariate_cols), with = FALSE],
  analysis_ids,
  "IID"
)
prs_dt <- align_to_ids(
  prs_raw_dt[, c("IID", prs_avg_cols), with = FALSE],
  analysis_ids,
  "IID"
)
exposure_dt <- align_to_ids(
  exposure_raw_dt[, c("IID", present_exposure_cols), with = FALSE],
  analysis_ids,
  "IID"
)

covariate_dt[, Age_z := zscore_strict(as.numeric(Age))]
covariate_dt[, avMSE_z := zscore_strict(as.numeric(avMSE))]
covariate_dt[, Sex := factor(Sex)]
covariate_dt[, gbatch := factor(gbatch)]

prs_z_dt <- copy(prs_dt)
for (prs_col in prs_avg_cols) {
  prs_z_dt[[prs_col]] <- zscore_strict(as.numeric(prs_z_dt[[prs_col]]))
}

analysis_id_master_dt <- data.table(IID = analysis_ids)
analysis_id_master_dt[, in_prs := IID %in% prs_raw_dt$IID]
analysis_id_master_dt[, in_oct_feature := IID %in% oct_feature_raw_dt$IID]
analysis_id_master_dt[, in_covariate := IID %in% covariate_raw_dt$IID]
analysis_id_master_dt[, in_exposure := IID %in% exposure_raw_dt$IID]
analysis_id_master_dt[, in_all_inputs := in_prs & in_oct_feature & in_covariate & in_exposure]

sample_alignment_summary_dt <- rbindlist(list(
  data.table(
    source = "prs",
    n_rows_raw = nrow(prs_raw_dt),
    n_unique_ids_raw = uniqueN(prs_raw_dt$IID),
    n_overlap_prs_ids = length(analysis_ids),
    n_aligned_rows = nrow(prs_dt),
    n_missing_rows_after_align = sum(!analysis_ids %in% prs_raw_dt$IID)
  ),
  data.table(
    source = "oct_feature",
    n_rows_raw = nrow(oct_feature_raw_dt),
    n_unique_ids_raw = uniqueN(oct_feature_raw_dt$IID),
    n_overlap_prs_ids = sum(analysis_ids %in% oct_feature_raw_dt$IID),
    n_aligned_rows = nrow(feature_dt),
    n_missing_rows_after_align = sum(!analysis_ids %in% oct_feature_raw_dt$IID)
  ),
  data.table(
    source = "covariate",
    n_rows_raw = nrow(covariate_raw_dt),
    n_unique_ids_raw = uniqueN(covariate_raw_dt$IID),
    n_overlap_prs_ids = sum(analysis_ids %in% covariate_raw_dt$IID),
    n_aligned_rows = nrow(covariate_dt),
    n_missing_rows_after_align = sum(!analysis_ids %in% covariate_raw_dt$IID)
  ),
  data.table(
    source = "exposure",
    n_rows_raw = nrow(exposure_raw_dt),
    n_unique_ids_raw = uniqueN(exposure_raw_dt$IID),
    n_overlap_prs_ids = sum(analysis_ids %in% exposure_raw_dt$IID),
    n_aligned_rows = nrow(exposure_dt),
    n_missing_rows_after_align = sum(!analysis_ids %in% exposure_raw_dt$IID)
  ),
  data.table(
    source = "full_intersection",
    n_rows_raw = NA_integer_,
    n_unique_ids_raw = NA_integer_,
    n_overlap_prs_ids = analysis_id_master_dt[in_all_inputs == TRUE, .N],
    n_aligned_rows = analysis_id_master_dt[in_all_inputs == TRUE, .N],
    n_missing_rows_after_align = 0L
  ),
  data.table(
    source = "step3_exposure_columns",
    n_rows_raw = nrow(step3_exposures),
    n_unique_ids_raw = NA_integer_,
    n_overlap_prs_ids = NA_integer_,
    n_aligned_rows = length(present_exposure_cols),
    n_missing_rows_after_align = length(missing_exposure_cols)
  )
), fill = TRUE)

sample_alignment_summary_dt

discovery_association_dt <- fread(discovery_association_path, check.names = FALSE)
sig_p_cutoff <- 0.05 / nrow(step3_exposures) / length(eye_features)

scan_grid_dt <- discovery_association_dt[
  term %in% step3_exposures$association_exposure &
    oct_features %in% eye_features,
  .(
    EyeFeature = oct_features,
    Exposure = term,
    discovery_beta = as.numeric(estimate),
    discovery_se = as.numeric(std.error),
    discovery_t = as.numeric(statistic),
    discovery_p = as.numeric(p.value),
    discovery_ci_low = as.numeric(conf.low),
    discovery_ci_high = as.numeric(conf.high),
    discovery_nobs = as.integer(nobs),
    discovery_p_bonferroni = as.numeric(p.value.bonferroni),
    discovery_p_fdr = as.numeric(p.value.fdr)
  )
]

scan_grid_dt <- merge(
  scan_grid_dt,
  step3_exposures,
  by.x = "Exposure",
  by.y = "association_exposure",
  all.x = TRUE,
  sort = FALSE
)
scan_grid_dt[, row_id := .I]
setcolorder(scan_grid_dt, "row_id")
scan_grid_dt[, is_primary := discovery_p < sig_p_cutoff]
stopifnot(nrow(scan_grid_dt) == 16704L, sum(scan_grid_dt$is_primary) == 640L)
scan_grid_dt[, primary_definition := primary_definition_value]
scan_grid_dt[, mapped_exposure_col := expression_exposure]
scan_grid_dt[, exposure_available_in_rep1 := mapped_exposure_col %in% present_exposure_cols]
scan_grid_dt[, prs_col := paste0(EyeFeature, "_AVG")]

processed_exposure_dt <- data.table(IID = exposure_dt$IID)
exposure_coding_list <- vector("list", length(present_exposure_cols))

for (exposure_idx in seq_along(present_exposure_cols)) {
  exposure_col <- present_exposure_cols[[exposure_idx]]
  exposure_type <- step3_exposures[expression_exposure == exposure_col, Type][1]
  semantic_numeric_opposite <- step3_exposures[
    expression_exposure == exposure_col,
    semantic_numeric_opposite
  ][1]
  processed <- process_exposure_vector(
    raw_vec = exposure_dt[[exposure_col]],
    exposure_name = exposure_col,
    exposure_type = exposure_type,
    semantic_numeric_opposite = semantic_numeric_opposite
  )
  processed_exposure_dt[[exposure_col]] <- processed$model_vec

  raw_num <- as.numeric(exposure_dt[[exposure_col]])
  raw_n_unique <- uniqueN(raw_num[!is.na(raw_num)])
  raw_min_group_n <- count_min_group_n(raw_num)

  exposure_coding_list[[exposure_idx]] <- data.table(
    expression_exposure = exposure_col,
    exposure_processed = processed$exposure_processed,
    exposure_original_coding = processed$exposure_original_coding,
    exposure_model_coding = processed$exposure_model_coding,
    semantic_transformation = processed$semantic_transformation,
    reference_level = processed$reference_level,
    exposed_level = processed$exposed_level,
    exposure_min = safe_min(raw_num),
    exposure_max = safe_max(raw_num),
    exposure_mean = safe_mean(raw_num),
    exposure_sd = safe_sd(raw_num),
    exposure_n_unique = raw_n_unique,
    exposure_level_counts = collapse_level_counts(raw_num),
    exposure_model_mean = safe_mean(processed$model_vec),
    exposure_model_sd = safe_sd(processed$model_vec),
    min_group_n = raw_min_group_n,
    sparse_exposure = exposure_type %in% c("Binary", "CatUnord") &&
      !is.na(raw_min_group_n) &&
      raw_min_group_n < 50
  )
}

missing_exposure_coding_dt <- step3_exposures[
  expression_exposure %in% missing_exposure_cols,
  .(
    expression_exposure,
    exposure_processed = "missing_exposure_in_rep1",
    exposure_original_coding = NA_character_,
    exposure_model_coding = NA_character_,
    semantic_transformation = NA_character_,
    reference_level = NA_character_,
    exposed_level = NA_character_,
    exposure_min = NA_real_,
    exposure_max = NA_real_,
    exposure_mean = NA_real_,
    exposure_sd = NA_real_,
    exposure_n_unique = NA_integer_,
    exposure_level_counts = NA_character_,
    exposure_model_mean = NA_real_,
    exposure_model_sd = NA_real_,
    min_group_n = NA_integer_,
    sparse_exposure = NA
  )
]

exposure_coding_dt <- rbindlist(
  c(exposure_coding_list, list(missing_exposure_coding_dt)),
  fill = TRUE
)

scan_grid_dt <- merge(
  scan_grid_dt,
  exposure_coding_dt,
  by = "expression_exposure",
  all.x = TRUE,
  sort = FALSE
)
setorder(scan_grid_dt, row_id)

scan_grid_dt[, .(
  n_rows = .N,
  n_primary = sum(is_primary),
  n_missing_exposure_in_rep1 = sum(!exposure_available_in_rep1)
)]

cov_model_cols <- c(
  "Age_z",
  "Sex",
  "avMSE_z",
  "gbatch",
  paste0("PC", 1:10)
)

fit_interaction_one <- function(scan_row) {
  scan_row <- as.data.table(scan_row)
  eye_feature <- scan_row$EyeFeature[[1]]
  exposure_col <- scan_row$expression_exposure[[1]]
  exposure_type <- scan_row$Type[[1]]
  prs_col <- scan_row$prs_col[[1]]

  base_result <- copy(scan_row)
  base_result[, `:=`(
    model_formula = paste0(
      eye_feature,
      " ~ Age_z + Sex + avMSE_z + gbatch + PC1 + ... + PC10 + ",
      prs_col,
      "_z + E_model + ",
      prs_col,
      "_z:E_model"
    ),
    n_complete = NA_integer_,
    n_missing_trait = NA_integer_,
    n_missing_prs = NA_integer_,
    n_missing_exposure = NA_integer_,
    n_missing_covariates = NA_integer_,
    beta_prs = NA_real_,
    se_prs = NA_real_,
    t_prs = NA_real_,
    p_prs = NA_real_,
    beta_exposure = NA_real_,
    se_exposure = NA_real_,
    t_exposure = NA_real_,
    p_exposure = NA_real_,
    beta_int = NA_real_,
    se_int = NA_real_,
    t_int = NA_real_,
    p_int = NA_real_,
    ci_low_int = NA_real_,
    ci_high_int = NA_real_,
    neglog10_p_int = NA_real_,
    direction_int = NA_character_,
    beta_prs_reference = NA_real_,
    beta_prs_exposed = NA_real_,
    prs_effect_E_minus1 = NA_real_,
    prs_effect_E_0 = NA_real_,
    prs_effect_E_plus1 = NA_real_,
    status = "ok",
    model_error = NA_character_
  )]

  if (!isTRUE(scan_row$exposure_available_in_rep1[[1]])) {
    base_result[, status := "missing_exposure_in_rep1"]
    return(base_result)
  }

  model_dt <- data.table(
    Y_trait = as.numeric(feature_dt[[eye_feature]]),
    Age_z = covariate_dt$Age_z,
    Sex = covariate_dt$Sex,
    avMSE_z = covariate_dt$avMSE_z,
    gbatch = covariate_dt$gbatch,
    PC1 = as.numeric(covariate_dt$PC1),
    PC2 = as.numeric(covariate_dt$PC2),
    PC3 = as.numeric(covariate_dt$PC3),
    PC4 = as.numeric(covariate_dt$PC4),
    PC5 = as.numeric(covariate_dt$PC5),
    PC6 = as.numeric(covariate_dt$PC6),
    PC7 = as.numeric(covariate_dt$PC7),
    PC8 = as.numeric(covariate_dt$PC8),
    PC9 = as.numeric(covariate_dt$PC9),
    PC10 = as.numeric(covariate_dt$PC10),
    PRS_trait_z = as.numeric(prs_z_dt[[prs_col]]),
    E_model = as.numeric(processed_exposure_dt[[exposure_col]])
  )

  cov_complete <- complete.cases(model_dt[, ..cov_model_cols])
  model_complete <- complete.cases(model_dt)

  base_result[, `:=`(
    n_missing_trait = sum(is.na(model_dt$Y_trait)),
    n_missing_prs = sum(is.na(model_dt$PRS_trait_z)),
    n_missing_exposure = sum(is.na(model_dt$E_model)),
    n_missing_covariates = sum(!cov_complete),
    n_complete = sum(model_complete)
  )]

  if (base_result$n_complete < 30) {
    base_result[, status := "insufficient_complete_cases"]
    return(base_result)
  }

  model_complete_dt <- model_dt[model_complete]

  if (uniqueN(model_complete_dt$E_model) < 2) {
    base_result[, status := "no_exposure_variation"]
    return(base_result)
  }

  if (uniqueN(model_complete_dt$PRS_trait_z) < 2) {
    base_result[, status := "no_prs_variation"]
    return(base_result)
  }

  fitted_model <- try(
    lm(
      Y_trait ~ Age_z + Sex + avMSE_z + gbatch +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
        PRS_trait_z * E_model,
      data = model_complete_dt
    ),
    silent = TRUE
  )

  if (inherits(fitted_model, "try-error")) {
    base_result[, `:=`(
      status = "model_error",
      model_error = as.character(fitted_model)
    )]
    return(base_result)
  }

  coef_dt <- as.data.table(coef(summary(fitted_model)), keep.rownames = "term")
  interaction_terms <- c("PRS_trait_z:E_model", "E_model:PRS_trait_z")
  interaction_term <- interaction_terms[interaction_terms %in% coef_dt$term][1]

  if (is.na(interaction_term)) {
    base_result[, status := "interaction_term_absent"]
    return(base_result)
  }

  beta_prs_value <- coef_dt[term == "PRS_trait_z", Estimate][1]
  se_prs_value <- coef_dt[term == "PRS_trait_z", `Std. Error`][1]
  t_prs_value <- coef_dt[term == "PRS_trait_z", `t value`][1]
  p_prs_value <- coef_dt[term == "PRS_trait_z", `Pr(>|t|)`][1]
  beta_exposure_value <- coef_dt[term == "E_model", Estimate][1]
  se_exposure_value <- coef_dt[term == "E_model", `Std. Error`][1]
  t_exposure_value <- coef_dt[term == "E_model", `t value`][1]
  p_exposure_value <- coef_dt[term == "E_model", `Pr(>|t|)`][1]
  beta_int_value <- coef_dt[term == interaction_term, Estimate][1]
  se_int_value <- coef_dt[term == interaction_term, `Std. Error`][1]
  t_int_value <- coef_dt[term == interaction_term, `t value`][1]
  p_int_value <- coef_dt[term == interaction_term, `Pr(>|t|)`][1]

  base_result[, `:=`(
    beta_prs = beta_prs_value,
    se_prs = se_prs_value,
    t_prs = t_prs_value,
    p_prs = p_prs_value,
    beta_exposure = beta_exposure_value,
    se_exposure = se_exposure_value,
    t_exposure = t_exposure_value,
    p_exposure = p_exposure_value,
    beta_int = beta_int_value,
    se_int = se_int_value,
    t_int = t_int_value,
    p_int = p_int_value,
    ci_low_int = beta_int_value - 1.96 * se_int_value,
    ci_high_int = beta_int_value + 1.96 * se_int_value,
    neglog10_p_int = -log10(p_int_value),
    direction_int = fifelse(beta_int_value > 0, "positive", "negative")
  )]

  if (exposure_type %in% c("Binary", "CatUnord")) {
    base_result[, `:=`(
      beta_prs_reference = beta_prs_value,
      beta_prs_exposed = beta_prs_value + beta_int_value
    )]
  } else {
    base_result[, `:=`(
      prs_effect_E_minus1 = beta_prs_value - beta_int_value,
      prs_effect_E_0 = beta_prs_value,
      prs_effect_E_plus1 = beta_prs_value + beta_int_value
    )]
  }

  base_result
}

fit_interaction_chunk <- function(chunk_dt) {
  rbindlist(
    lapply(seq_len(nrow(chunk_dt)), function(i) fit_interaction_one(chunk_dt[i])),
    fill = TRUE
  )
}

scan_chunks <- split(scan_grid_dt, by = "EyeFeature", keep.by = TRUE)
interaction_result_dt <- rbindlist(
  mclapply(
    scan_chunks,
    fit_interaction_chunk,
    mc.cores = min(n_cores, length(scan_chunks))
  ),
  fill = TRUE
)
setorder(interaction_result_dt, row_id)

main_effect_filter_pass_idx <- interaction_result_dt$status == "ok" &
  !is.na(interaction_result_dt$p_prs) &
  !is.na(interaction_result_dt$p_exposure) &
  !is.na(interaction_result_dt$p_int) &
  interaction_result_dt$p_prs < 0.05 &
  interaction_result_dt$p_exposure < 0.05

main_effect_filter_all_idx <- main_effect_filter_pass_idx
main_effect_filter_primary_idx <- main_effect_filter_pass_idx &
  interaction_result_dt$is_primary == TRUE

interaction_result_dt[, `:=`(
  main_effect_filter_pass = main_effect_filter_pass_idx,
  p_int_fdr_main_effect_filter_all = NA_real_,
  p_int_bonf_main_effect_filter_all = NA_real_,
  p_int_fdr_main_effect_filter_primary = NA_real_,
  p_int_bonf_main_effect_filter_primary = NA_real_
)]

if (sum(main_effect_filter_all_idx) > 0) {
  interaction_result_dt[main_effect_filter_all_idx, `:=`(
    p_int_fdr_main_effect_filter_all = p.adjust(p_int, method = "BH"),
    p_int_bonf_main_effect_filter_all = pmin(p_int * .N, 1)
  )]
}

if (sum(main_effect_filter_primary_idx) > 0) {
  interaction_result_dt[main_effect_filter_primary_idx, `:=`(
    p_int_fdr_main_effect_filter_primary = p.adjust(p_int, method = "BH"),
    p_int_bonf_main_effect_filter_primary = pmin(p_int * .N, 1)
  )]
}

interaction_result_dt[, `:=`(
  significance_fdr_main_effect_filter_all = !is.na(p_int_fdr_main_effect_filter_all) &
    p_int_fdr_main_effect_filter_all < 0.05,
  significance_bonf_main_effect_filter_all = !is.na(p_int_bonf_main_effect_filter_all) &
    p_int_bonf_main_effect_filter_all < 0.05,
  significance_fdr_main_effect_filter_primary = !is.na(p_int_fdr_main_effect_filter_primary) &
    p_int_fdr_main_effect_filter_primary < 0.05,
  significance_bonf_main_effect_filter_primary = !is.na(p_int_bonf_main_effect_filter_primary) &
    p_int_bonf_main_effect_filter_primary < 0.05,
  analysis_type = fifelse(is_primary, "primary", "exploratory")
)]

interaction_annotated_dt <- merge(
  interaction_result_dt,
  oct_dictionary_dt,
  by = "EyeFeature",
  all.x = TRUE,
  sort = FALSE
)
setorder(interaction_annotated_dt, row_id)

main_effect_filtered_dt <- interaction_annotated_dt[main_effect_filter_pass == TRUE]
setorder(main_effect_filtered_dt, p_int)

top_filtered_interactions_dt <- head(main_effect_filtered_dt, 100)

summary_dt <- rbindlist(list(
  data.table(metric = "n_scan_rows", value = nrow(interaction_annotated_dt)),
  data.table(metric = "n_status_ok", value = interaction_annotated_dt[status == "ok", .N]),
  data.table(metric = "n_primary_rows", value = interaction_annotated_dt[is_primary == TRUE, .N]),
  data.table(metric = "n_primary_status_ok", value = interaction_annotated_dt[is_primary == TRUE & status == "ok", .N]),
  data.table(metric = "n_main_effect_filter_pass", value = interaction_annotated_dt[main_effect_filter_pass == TRUE, .N]),
  data.table(metric = "n_primary_main_effect_filter_pass", value = interaction_annotated_dt[main_effect_filter_pass == TRUE & is_primary == TRUE, .N]),
  data.table(metric = "n_fdr_main_effect_filter_all_0.05", value = interaction_annotated_dt[significance_fdr_main_effect_filter_all == TRUE, .N]),
  data.table(metric = "n_bonf_main_effect_filter_all_0.05", value = interaction_annotated_dt[significance_bonf_main_effect_filter_all == TRUE, .N]),
  data.table(metric = "n_fdr_main_effect_filter_primary_0.05", value = interaction_annotated_dt[significance_fdr_main_effect_filter_primary == TRUE, .N]),
  data.table(metric = "n_bonf_main_effect_filter_primary_0.05", value = interaction_annotated_dt[significance_bonf_main_effect_filter_primary == TRUE, .N]),
  data.table(metric = "n_sparse_exposure_rows", value = interaction_annotated_dt[sparse_exposure == TRUE, .N]),
  data.table(metric = "n_missing_exposure_in_rep1_rows", value = interaction_annotated_dt[status == "missing_exposure_in_rep1", .N]),
  interaction_annotated_dt[, .(value = .N), by = .(metric = paste0("status__", status))],
  interaction_annotated_dt[, .(value = .N), by = .(metric = paste0("type__", Type))],
  interaction_annotated_dt[, .(value = .N), by = .(metric = paste0("category__", gsub("[^A-Za-z0-9]+", "_", Category)))]
), fill = TRUE)

write_tsv(interaction_annotated_dt, interaction_with_main_effect_p_path)
write_tsv(main_effect_filtered_dt, main_effect_filtered_path)
write_tsv(summary_dt, main_effect_filter_summary_path)
write_tsv(top_filtered_interactions_dt, top_filtered_interactions_path)

summary_dt
