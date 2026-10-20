# Fixed selected tasks; no concordant-pair or exposure selection is performed here.
library(data.table)
library(survival)
data.table::setDTthreads(2L)
task_table <- data.table::fread("data/selected_disease_exposure_tasks.csv")
prs_matrix_path <- "data/white_non_discovery_prs.tsv"
covariate_path <- "data/disease_covariates.csv"
pca_path <- "data/disease_pca.txt"
exposure_matrix_path <- "data/raw_disease_exposures.csv"
disease_dir <- "data/incident_diseases"
result_dir <- "results/prs_stratified_disease"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
result_path <- file.path(result_dir, "stratified_cox_results.csv")
group_summary_path <- file.path(result_dir, "group_summary.csv")
prs_cutpoint_path <- file.path(result_dir, "prs_quartile_cutpoints.csv")
exposure_cutpoint_path <- file.path(result_dir, "exposure_cutpoints.csv")
task_manifest_path <- file.path(result_dir, "task_manifest.csv")
qc_path <- file.path(result_dir, "model_status.csv")
expected_task_count <- 109L
max_workers <- 50L
runtime_workers <- 50L
prs_group_levels <- paste0("Q", 1:4)
pc_columns <- paste0("PC", 1:10)
covariate_terms <- paste(c("Age", "Sex", "avMSE_imputed", "gbatch", pc_columns), collapse = " + ")
# 读取PRS、协变量和PC；PRS先在完整PRS人群中标准化。
prs_columns <- paste0(unique(task_table$retinal_feature), "_AVG")
prs_data <- fread(prs_matrix_path, select = c("IID", prs_columns))
prs_data[, IID := as.character(IID)]
if (anyDuplicated(prs_data$IID)) stop("Duplicated PRS IID")
prs_data[, (prs_columns) := lapply(.SD, function(x) as.numeric(scale(as.numeric(x)))), .SDcols = prs_columns]

covariates <- fread(covariate_path, select = c("eid", "Age", "Sex", "avMSE", "gbatch"))
covariates[, eid := as.character(eid)]
if (anyDuplicated(covariates$eid)) stop("Duplicated covariate eid")
pca_data <- as.data.table(utils::read.table(
  pca_path,
  header = FALSE,
  sep = "",
  fill = TRUE,
  col.names = c("eid", paste0("PC", 1:40)),
  colClasses = c("character", rep("numeric", 10), rep("NULL", 30))
))
pca_data[, eid := as.character(eid)]
if (anyDuplicated(pca_data$eid)) stop("Duplicated PCA eid")

aligned_covariates <- covariates[match(prs_data$IID, covariates$eid), .(Age, Sex, avMSE, gbatch)]
aligned_pca <- pca_data[match(prs_data$IID, pca_data$eid), .SD, .SDcols = pc_columns]
prs_data[, c("Age", "Sex", "avMSE", "gbatch") := aligned_covariates]
prs_data[, (pc_columns) := aligned_pca]

# 读取全人群暴露矩阵；只保留109个任务实际使用的表达列。
exposure_ids <- unique(task_table$expression_exposure)
exposure_data <- fread(exposure_matrix_path, select = c("eid", exposure_ids))
exposure_data[, eid := as.character(eid)]
if (anyDuplicated(exposure_data$eid)) stop("Duplicated exposure eid")

# 读取疾病结局。
disease_levels <- unique(task_table$disease)
disease_data <- setNames(lapply(disease_levels, function(disease_name) {
  x <- fread(file.path(disease_dir, paste0(disease_name, ".csv")), select = c("eid", "BL2Target_yrs", "target_y"))
  x[, eid := as.character(eid)]
  if (anyDuplicated(x$eid)) stop(paste("Duplicated disease eid:", disease_name))
  x
}), disease_levels)

# 为每个疾病–表型组合建立共同基础队列，并固定PRS四分位、avMSE填补和暴露三分位口径。
pair_ids <- paste(task_table$disease, task_table$retinal_feature, sep = "__")
task_table[, pair_id := paste(disease, retinal_feature, sep = "__")]
pair_bases <- list()
prs_cutpoints <- list()
exposure_cutpoints <- list()

for (pair_id in unique(task_table$pair_id)) {
  pair_value <- pair_id
  pair_task <- task_table[pair_id == pair_value]
  disease_name <- pair_task$disease[1L]
  retinal_feature <- pair_task$retinal_feature[1L]
  prs_id <- paste0(retinal_feature, "_AVG")
  outcome <- disease_data[[disease_name]]
  outcome_index <- match(prs_data$IID, outcome$eid)
  base <- data.table(
    eid = prs_data$IID,
    PRS_z = prs_data[[prs_id]],
    Age = as.numeric(prs_data$Age),
    Sex = factor(prs_data$Sex),
    avMSE = as.numeric(prs_data$avMSE),
    gbatch = factor(prs_data$gbatch)
  )
  base[, (pc_columns) := prs_data[, .SD, .SDcols = pc_columns]]
  base[, time := outcome$BL2Target_yrs[outcome_index]]
  base[, event := outcome$target_y[outcome_index]]
  base <- base[
    is.finite(PRS_z) & !is.na(Age) & !is.na(Sex) & !is.na(gbatch) &
      complete.cases(base[, .SD, .SDcols = pc_columns]) & !is.na(time) & !is.na(event)
  ]
  if (nrow(base) == 0L) stop(paste("Empty base cohort:", pair_id))
  base[, avMSE_was_missing := is.na(avMSE)]
  base[, avMSE_imputed := {
    sex_median <- median(avMSE, na.rm = TRUE)
    replace(avMSE, is.na(avMSE), sex_median)
  }, by = Sex]
  if (anyNA(base$avMSE_imputed)) stop(paste("avMSE imputation left missing values:", pair_id))
  prs_cuts <- quantile(base$PRS_z, probs = c(.25, .5, .75), names = FALSE, na.rm = TRUE)
  if (any(!is.finite(prs_cuts)) || anyDuplicated(prs_cuts)) stop(paste("Invalid PRS quartile cutpoints:", pair_id))
  base[, prs_group := cut(PRS_z, breaks = c(-Inf, prs_cuts, Inf), labels = prs_group_levels, include.lowest = TRUE)]
  pair_exposure_ids <- unique(pair_task$expression_exposure)
  base[, (pair_exposure_ids) := exposure_data[match(base$eid, exposure_data$eid), .SD, .SDcols = pair_exposure_ids]]
  pair_bases[[pair_id]] <- base
  prs_cutpoints[[pair_id]] <- data.table(
    pair_id, disease = disease_name, retinal_feature,
    boundary = c("Q2/Q1", "Q3/Q2", "Q4/Q3"), cutpoint = prs_cuts,
    n_base = nrow(base), n_event_base = sum(base$event == 1L)
  )
  exposure_cutpoints[[pair_id]] <- rbindlist(lapply(pair_task$expression_exposure, function(exposure_id) {
    map_row <- pair_task[expression_exposure == exposure_id][1L]
    if (map_row$exposure_type %in% c("CatOrd", "Cont")) {
      raw <- as.numeric(base[[exposure_id]])
      if (isTRUE(as.logical(map_row$semantic_numeric_opposite))) raw <- 7 - raw
      raw <- raw[is.finite(raw)]
      if (length(raw) == 0L) stop(paste("Invalid exposure tertile cutpoints:", pair_id, exposure_id))
      ordered_values <- sort(unique(raw))
      if (length(ordered_values) < 2L) stop(paste("Exposure has fewer than two levels:", pair_id, exposure_id))
      cuts <- quantile(raw, probs = c(1 / 3, 2 / 3), names = FALSE, na.rm = TRUE)
      tertile_method <- "quantile_tertile"
      if (length(ordered_values) == 2L) {
        cuts <- ordered_values[1L]
        tertile_method <- "ordered_levels_two_group_fallback"
      } else if (any(!is.finite(cuts)) || anyDuplicated(cuts)) {
        level_counts <- tabulate(match(raw, ordered_values), nbins = length(ordered_values))
        cumulative_fraction <- cumsum(level_counts) / sum(level_counts)
        boundary_pairs <- expand.grid(
          first = seq_len(length(ordered_values) - 1L),
          second = seq_len(length(ordered_values) - 1L)
        )
        boundary_pairs <- boundary_pairs[boundary_pairs$first < boundary_pairs$second, , drop = FALSE]
        boundary_pairs$loss <-
          (cumulative_fraction[boundary_pairs$first] - 1 / 3)^2 +
          (cumulative_fraction[boundary_pairs$second] - 2 / 3)^2
        best_pair <- boundary_pairs[which.min(boundary_pairs$loss), c("first", "second")]
        cuts <- ordered_values[c(best_pair$first, best_pair$second)]
        tertile_method <- "ordered_levels_balanced_tertile"
      }
      # 切点等于最大值时High为空：删除空区间的边界，保留既有非空组。
      group_counts <- tabulate(cut(raw, c(-Inf, cuts, Inf), labels = FALSE), nbins = length(cuts) + 1L)
      if (length(cuts) == 2L && any(group_counts == 0L)) {
        cuts <- cuts[vapply(cuts, function(x) any(raw <= x) && any(raw > x), logical(1))]
        if (length(cuts) == 2L && any(group_counts == 0L)) cuts <- cuts[1L]
        if (length(cuts) != 1L) stop(paste("Cannot form two nonempty groups:", pair_id, exposure_id))
        tertile_method <- "ordered_levels_two_group_fallback"
      }
      if (any(!is.finite(cuts)) || anyDuplicated(cuts)) stop(paste("Invalid exposure tertile cutpoints:", pair_id, exposure_id))
      data.table(pair_id, disease = disease_name, retinal_feature, expression_exposure = exposure_id,
        boundary = if (length(cuts) == 1L) "High/Low" else c("Medium/Low", "High/Medium"), cutpoint = cuts,
        exposure_transform = fifelse(isTRUE(as.logical(map_row$semantic_numeric_opposite)), "7_minus_x", "none"),
        tertile_method)
    } else {
      data.table(pair_id, disease = disease_name, retinal_feature, expression_exposure = exposure_id,
        boundary = NA_character_, cutpoint = NA_real_, exposure_transform = "none", tertile_method = "not_applicable")
    }
  }), fill = TRUE)
}
names(pair_bases) <- unique(task_table$pair_id)
prs_cutpoints <- rbindlist(prs_cutpoints, fill = TRUE)
exposure_cutpoints <- rbindlist(exposure_cutpoints, fill = TRUE)

# 单个PRS四分位和单个暴露任务的Cox拟合。
fit_task <- function(task_index) {
  task_row <- task_table[task_index]
  pair_base <- data.table::copy(pair_bases[[task_row$pair_id]])
  exposure_id <- task_row$expression_exposure
  raw <- as.numeric(pair_base[[exposure_id]])
  exposure_transform <- "none"
  if (isTRUE(as.logical(task_row$semantic_numeric_opposite))) {
    raw <- 7 - raw
    exposure_transform <- "7_minus_x"
  }
  if (task_row$exposure_type %in% c("Binary", "CatUnord")) {
    if (any(!is.na(raw) & !raw %in% c(0, 1))) stop(paste("Non-binary exposure values:", exposure_id))
    exposure_group_local <- fifelse(is.na(raw), NA_character_, fifelse(raw == 0, "reference", "exposed"))
    exposure_levels <- c("reference", "exposed")
    reference_level <- "reference"
    display_labels <- c(reference = "Reference", exposed = "Exposed")
  } else {
    cuts <- exposure_cutpoints[pair_id == task_row$pair_id & expression_exposure == exposure_id & !is.na(cutpoint), cutpoint]
    exposure_levels <- if (length(cuts) == 1L) c("low", "high") else c("low", "medium", "high")
    exposure_group_local <- cut(raw, breaks = c(-Inf, cuts, Inf), labels = exposure_levels, include.lowest = TRUE)
    reference_level <- "low"
    display_labels <- c(low = "Low", medium = "Medium", high = "High")
  }
  pair_base[, exposure_group := factor(exposure_group_local, levels = exposure_levels)]
  expected_contrasts <- setdiff(exposure_levels, reference_level)
  missing_exposure <- sum(is.na(pair_base$exposure_group))
  fit_results <- list()
  group_results <- list()
  failure_messages <- character()
  all_warning_messages <- character()
  for (prs_level in prs_group_levels) {
    model_data <- pair_base[prs_group == prs_level & !is.na(exposure_group)]
    model_data[, exposure_group := droplevels(exposure_group)]
    group_results[[prs_level]] <- model_data[, .(
      disease = task_row$disease, retinal_feature = task_row$retinal_feature,
      expression_exposure = exposure_id, prs_group = prs_level,
      n = .N, n_event = sum(event == 1L)
    ), by = exposure_group]
    warning_messages <- character()
    fit_error <- NULL
    fit <- tryCatch(
      withCallingHandlers(
        coxph(as.formula(paste("Surv(time, event) ~ exposure_group +", covariate_terms)),
          data = model_data, ties = "efron", na.action = na.fail),
        warning = function(w) {
          warning_messages <<- c(warning_messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
    all_warning_messages <- c(all_warning_messages, warning_messages)
    coefficients <- if (!is.null(fit)) coef(fit) else numeric()
    covariance <- if (!is.null(fit)) vcov(fit) else matrix(numeric(), 0, 0)
    for (contrast in expected_contrasts) {
      term <- paste0("exposure_group", contrast)
      status <- "fitted"
      error_text <- NA_character_
      if (!is.null(fit_error)) {
        status <- "failed"
        error_text <- fit_error
      } else if (!(term %in% names(coefficients)) || !is.finite(coefficients[[term]]) || !is.finite(covariance[term, term])) {
        status <- "failed"
        error_text <- paste("Missing or non-finite coefficient:", term)
      }
      if (status == "fitted") {
        beta <- coefficients[[term]]
        se <- sqrt(covariance[term, term])
        fit_results[[length(fit_results) + 1L]] <- data.table(
          model_status = status, task_id = task_row$task_id, pair_id = task_row$pair_id,
          disease = task_row$disease, retinal_feature = task_row$retinal_feature,
          feature_name = task_row$feature_name, exposure_status = task_row$exposure_status,
          association_exposure = task_row$association_exposure, expression_exposure = exposure_id,
          exposure_name = task_row$exposure_name, exposure_category = task_row$exposure_category,
          exposure_type = task_row$exposure_type, exposure_transform,
          prs_group = prs_level, exposure_group = display_labels[[contrast]], reference_group = display_labels[[reference_level]],
          beta, hazard_ratio = exp(beta), ci_low = exp(beta - qnorm(.975) * se), ci_high = exp(beta + qnorm(.975) * se),
          p_value = 2 * pnorm(-abs(beta / se)), n = nrow(model_data), n_event = sum(model_data$event == 1L),
          n_missing_exposure = missing_exposure, n_avMSE_imputed = sum(model_data$avMSE_was_missing),
          warning_text = paste(unique(warning_messages), collapse = " | "), error_text = NA_character_
        )
      } else {
        failure_messages <- c(failure_messages, paste(prs_level, contrast, error_text))
        fit_results[[length(fit_results) + 1L]] <- data.table(
          model_status = status, task_id = task_row$task_id, pair_id = task_row$pair_id,
          disease = task_row$disease, retinal_feature = task_row$retinal_feature,
          feature_name = task_row$feature_name, exposure_status = task_row$exposure_status,
          association_exposure = task_row$association_exposure, expression_exposure = exposure_id,
          exposure_name = task_row$exposure_name, exposure_category = task_row$exposure_category,
          exposure_type = task_row$exposure_type, exposure_transform,
          prs_group = prs_level, exposure_group = display_labels[[contrast]], reference_group = display_labels[[reference_level]],
          beta = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
          n = nrow(model_data), n_event = sum(model_data$event == 1L), n_missing_exposure = missing_exposure,
          n_avMSE_imputed = sum(model_data$avMSE_was_missing), warning_text = paste(unique(warning_messages), collapse = " | "),
          error_text
        )
      }
    }
  }
  manifest <- data.table(
    task_id = task_row$task_id, pair_id = task_row$pair_id, disease = task_row$disease,
    retinal_feature = task_row$retinal_feature, feature_name = task_row$feature_name,
    association_exposure = task_row$association_exposure, expression_exposure = exposure_id,
    exposure_name = task_row$exposure_name, exposure_status = task_row$exposure_status,
    exposure_type = task_row$exposure_type, exposure_transform,
    n_base = nrow(pair_base), n_event_base = sum(pair_base$event == 1L),
    n_missing_exposure = missing_exposure, n_avMSE_imputed_base = sum(pair_base$avMSE_was_missing),
    model_status = fifelse(length(failure_messages) == 0L, "fitted", "failed"),
    warning_text = if (length(all_warning_messages) == 0L) NA_character_ else paste(unique(all_warning_messages), collapse = " | "),
    error_text = if (length(failure_messages) == 0L) NA_character_ else paste(failure_messages, collapse = " | ")
  )
  list(
    results = rbindlist(fit_results, fill = TRUE),
    group_summary = rbindlist(group_results, fill = TRUE),
    manifest = manifest
  )
}

# 单个任务发生输入或模型错误时保留任务级失败记录，不把失败当作阴性结果。
fit_task_safe <- function(task_index) {
  tryCatch(
    fit_task(task_index),
    error = function(e) {
      task_row <- task_table[task_index]
      list(
        results = data.table(
          model_status = "failed", task_id = task_row$task_id, pair_id = task_row$pair_id,
          disease = task_row$disease, retinal_feature = task_row$retinal_feature,
          feature_name = task_row$feature_name, exposure_status = task_row$exposure_status,
          association_exposure = task_row$association_exposure, expression_exposure = task_row$expression_exposure,
          exposure_name = task_row$exposure_name, exposure_category = task_row$exposure_category,
          exposure_type = task_row$exposure_type, exposure_transform = NA_character_,
          prs_group = NA_character_, exposure_group = NA_character_, reference_group = NA_character_,
          beta = NA_real_, hazard_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
          n = NA_integer_, n_event = NA_integer_, n_missing_exposure = NA_integer_, n_avMSE_imputed = NA_integer_,
          warning_text = NA_character_, error_text = conditionMessage(e)
        ),
        group_summary = data.table(),
        manifest = data.table(
          task_id = task_row$task_id, pair_id = task_row$pair_id, disease = task_row$disease,
          retinal_feature = task_row$retinal_feature, feature_name = task_row$feature_name,
          association_exposure = task_row$association_exposure, expression_exposure = task_row$expression_exposure,
          exposure_name = task_row$exposure_name, exposure_status = task_row$exposure_status,
          exposure_type = task_row$exposure_type, exposure_transform = NA_character_,
          n_base = NA_integer_, n_event_base = NA_integer_, n_missing_exposure = NA_integer_,
          n_avMSE_imputed_base = NA_integer_, model_status = "failed",
          warning_text = NA_character_, error_text = conditionMessage(e)
        )
      )
    }
  )
}

# 按批次并行，保留进度和预计剩余时间。
fit_all <- vector("list", nrow(task_table))
analysis_start <- Sys.time()
for (batch_start in seq(1L, nrow(task_table), by = max_workers)) {
  batch_index <- batch_start:min(batch_start + max_workers - 1L, nrow(task_table))
  fit_all[batch_index] <- parallel::mclapply(
    batch_index,
    fit_task_safe,
    mc.cores = min(runtime_workers, length(batch_index)),
    mc.preschedule = FALSE
  )
  completed <- max(batch_index)
  elapsed_seconds <- as.numeric(difftime(Sys.time(), analysis_start, units = "secs"))
  eta_seconds <- if (completed > 0L) elapsed_seconds / completed * (nrow(task_table) - completed) else NA_real_
  cat(sprintf(
    "completed=%d/%d (%.1f%%); elapsed=%.1fs; ETA=%.1fs\n",
    completed, nrow(task_table), 100 * completed / nrow(task_table), elapsed_seconds, eta_seconds
  ))
}

stratified_results <- rbindlist(lapply(fit_all, `[[`, "results"), fill = TRUE)
group_summary <- rbindlist(lapply(fit_all, `[[`, "group_summary"), fill = TRUE)
task_manifest <- rbindlist(lapply(fit_all, `[[`, "manifest"), fill = TRUE)
if (nrow(task_manifest) != expected_task_count) stop("Task manifest does not contain 109 tasks")
if (any(stratified_results$model_status == "fitted" & !is.finite(stratified_results$hazard_ratio))) stop("Non-finite fitted hazard ratio")

qc_results <- task_manifest[, .(
  task_id, pair_id, disease, retinal_feature, expression_exposure, exposure_status,
  model_status, n_base, n_event_base, n_missing_exposure, n_avMSE_imputed_base,
  warning_text, error_text
)]

fwrite(stratified_results, result_path)
fwrite(group_summary, group_summary_path)
fwrite(prs_cutpoints, prs_cutpoint_path)
fwrite(exposure_cutpoints, exposure_cutpoint_path)
fwrite(task_manifest, task_manifest_path)
fwrite(qc_results, qc_path)
