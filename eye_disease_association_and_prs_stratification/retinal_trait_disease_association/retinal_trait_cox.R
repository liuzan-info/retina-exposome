# Discovery retinal-trait Cox associations, using precomputed incident endpoints.
options(stringsAsFactors = FALSE)
library(survival)
discovery_eid_path <- "data/discovery_eids.txt"
disease_dir <- "data/incident_diseases"
retinal_feature_matrix_path <- "data/retinal_phenotypes.tsv"
retinal_feature_dictionary_path <- "data/retinal_dictionary.csv"
covariate_path <- "data/gwas_covariates.tsv"
output_dir <- "results/retinal_trait_disease"
result_path <- file.path(output_dir, "cox_retinal_feature_eye_disease_results.csv")
sample_summary_path <- file.path(output_dir, "sample_summary.csv")
disease_names <- c("AMD", "RD", "Cataract", "POAG", "DR")
n_threads <- 10L
expected_incident_events <- c(
  AMD = 682L,
  RD = 270L,
  Cataract = 3463L,
  POAG = 659L,
  DR = 259L
)

required_disease_columns <- c(
  "eid", "target_y", "BL2Target_yrs", "case_timing"
)
required_dictionary_columns <- c(
  "EyeFeature",
  "feature_name",
  "plot_label_med",
  "group",
  "rank",
  "group_rank",
  "bubble",
  "group_color"
)

# 校验必需字段
check_required_columns <- function(data, required_columns, data_name) {
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0L) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

# 读取 discovery EID
discovery_eid <- scan(
  discovery_eid_path,
  what = character(),
  quiet = TRUE
)
if (length(discovery_eid) != 33030L) {
  stop("Discovery EID count is not 33,030: ", length(discovery_eid))
}
if (anyDuplicated(discovery_eid) > 0L) {
  stop("Discovery EID file contains duplicated EID")
}

# 读取 58 个眼底表型词典
retinal_feature_dictionary <- read.csv(
  retinal_feature_dictionary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)
colnames(retinal_feature_dictionary)[1] <- sub(
  "^\\ufeff",
  "",
  colnames(retinal_feature_dictionary)[1]
)
check_required_columns(
  retinal_feature_dictionary,
  required_dictionary_columns,
  "Retinal feature dictionary"
)
if (nrow(retinal_feature_dictionary) != 58L) {
  stop("Retinal feature dictionary does not contain 58 rows")
}
if (anyDuplicated(retinal_feature_dictionary$EyeFeature) > 0L) {
  stop("Retinal feature dictionary contains duplicated EyeFeature")
}
retinal_feature_names <- retinal_feature_dictionary$EyeFeature

# 读取上游标准化眼底表型矩阵，不再次标准化
retinal_feature_matrix <- read.table(
  retinal_feature_matrix_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
check_required_columns(
  retinal_feature_matrix,
  c("IID", retinal_feature_names),
  "Retinal feature matrix"
)
retinal_feature_matrix$IID <- as.character(retinal_feature_matrix$IID)
if (anyDuplicated(retinal_feature_matrix$IID) > 0L) {
  stop("Retinal feature matrix contains duplicated IID")
}
if (any(!discovery_eid %in% retinal_feature_matrix$IID)) {
  stop("Some discovery EID are missing from retinal feature matrix")
}
non_numeric_features <- retinal_feature_names[
  !vapply(
    retinal_feature_matrix[retinal_feature_names],
    is.numeric,
    logical(1)
  )
]
if (length(non_numeric_features) > 0L) {
  stop(
    "Retinal feature columns are not numeric: ",
    paste(non_numeric_features, collapse = ", ")
  )
}

# 读取协变量
covariate_data <- read.table(
  covariate_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
check_required_columns(
  covariate_data,
  c("IID", "Age", "Sex", "avMSE"),
  "Covariate data"
)
covariate_data$IID <- as.character(covariate_data$IID)
if (anyDuplicated(covariate_data$IID) > 0L) {
  stop("Covariate data contains duplicated IID")
}

# 读取单个 incident 疾病文件并限定 discovery 样本
read_discovery_disease <- function(disease_name) {
  disease_path <- file.path(disease_dir, paste0(disease_name, ".csv"))
  disease_data <- read.csv(
    disease_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  check_required_columns(
    disease_data,
    required_disease_columns,
    paste0(disease_name, " disease data")
  )

  disease_data$eid <- as.character(disease_data$eid)
  disease_data$target_y <- as.numeric(disease_data$target_y)
  disease_data$BL2Target_yrs <- as.numeric(disease_data$BL2Target_yrs)

  if (anyDuplicated(disease_data$eid) > 0L) {
    stop(disease_name, " incident data contains duplicated EID")
  }
  if (any(is.na(disease_data$target_y)) ||
      !all(disease_data$target_y %in% c(0, 1))) {
    stop(disease_name, " target_y contains missing or non-0/1 values")
  }
  if (any(is.na(disease_data$BL2Target_yrs)) ||
      any(disease_data$BL2Target_yrs < 0)) {
    stop(disease_name, " incident data contains missing or negative follow-up")
  }
  if (any(disease_data$case_timing == "prevalent")) {
    stop(disease_name, " incident data still contains prevalent cases")
  }

  disease_data <- disease_data[
    disease_data$eid %in% discovery_eid,
    ,
    drop = FALSE
  ]
  rownames(disease_data) <- NULL

  observed_events <- sum(disease_data$target_y == 1L)
  if (observed_events != expected_incident_events[[disease_name]]) {
    stop(
      disease_name,
      " discovery incident event count is ",
      observed_events,
      "; expected ",
      expected_incident_events[[disease_name]]
    )
  }

  disease_data
}

disease_data_list <- lapply(disease_names, read_discovery_disease)
names(disease_data_list) <- disease_names

model_tasks <- expand.grid(
  disease_name = disease_names,
  feature_index = seq_along(retinal_feature_names),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# 拟合单个 disease × retinal feature 模型
fit_retinal_feature_cox <- function(task_index) {
  task <- model_tasks[task_index, , drop = FALSE]
  disease_name <- task$disease_name[[1]]
  feature_index_in_dictionary <- task$feature_index[[1]]
  disease_data <- disease_data_list[[disease_name]]
  covariate_index <- match(disease_data$eid, covariate_data$IID)
  feature_index <- match(disease_data$eid, retinal_feature_matrix$IID)

  if (anyNA(covariate_index)) {
    stop(disease_name, " discovery disease EID are missing from covariate data")
  }
  if (anyNA(feature_index)) {
    stop(disease_name, " discovery disease EID are missing from retinal feature matrix")
  }

  retinal_feature_name <- retinal_feature_names[[feature_index_in_dictionary]]
  feature_metadata <- retinal_feature_dictionary[
    feature_index_in_dictionary,
    ,
    drop = FALSE
  ]

    model_data <- data.frame(
      BL2Target_yrs = disease_data$BL2Target_yrs,
      target_y = disease_data$target_y,
      predictor = retinal_feature_matrix[[retinal_feature_name]][feature_index],
      Age = covariate_data$Age[covariate_index],
      Sex = covariate_data$Sex[covariate_index],
      avMSE = covariate_data$avMSE[covariate_index],
      stringsAsFactors = FALSE
    )
    model_data <- model_data[
      complete.cases(model_data),
      ,
      drop = FALSE
    ]
    model_data$Sex <- droplevels(factor(model_data$Sex))

    if (length(unique(model_data$predictor)) < 2L) {
      stop(
        disease_name,
        " / ",
        retinal_feature_name,
        ": predictor has no variation"
      )
    }
    if (length(unique(model_data$target_y)) < 2L) {
      stop(
        disease_name,
        " / ",
        retinal_feature_name,
        ": outcome has no variation"
      )
    }
    if (nlevels(model_data$Sex) < 2L) {
      stop(
        disease_name,
        " / ",
        retinal_feature_name,
        ": Sex has fewer than two levels"
      )
    }

    cox_model <- survival::coxph(
      survival::Surv(BL2Target_yrs, target_y) ~
        predictor + Age + Sex + avMSE,
      data = model_data,
      ties = "efron",
      na.action = na.fail,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    coefficient_table <- summary(cox_model)$coefficients

    beta <- unname(coefficient_table["predictor", "coef"])
    standard_error <- unname(coefficient_table["predictor", "se(coef)"])
    z_value <- unname(coefficient_table["predictor", "z"])
    p_value <- unname(coefficient_table["predictor", "Pr(>|z|)"])

    if (!all(is.finite(c(beta, standard_error, z_value, p_value)))) {
      stop(
        disease_name,
        " / ",
        retinal_feature_name,
        ": non-finite Cox estimate"
      )
    }

    data.frame(
      disease = disease_name,
      retinal_feature = retinal_feature_name,
      feature_name = feature_metadata$feature_name,
      plot_label = feature_metadata$plot_label_med,
      feature_group = feature_metadata$group,
      feature_rank = feature_metadata$rank,
      group_rank = feature_metadata$group_rank,
      bubble = feature_metadata$bubble,
      group_color = feature_metadata$group_color,
      effect_unit = "per 1 upstream-normalized unit",
      covariates = "Age + Sex + avMSE",
      beta = beta,
      standard_error = standard_error,
      hazard_ratio = exp(beta),
      confidence_interval_low = exp(beta - qnorm(0.975) * standard_error),
      confidence_interval_high = exp(beta + qnorm(0.975) * standard_error),
      z_value = z_value,
      p_value = p_value,
      n = nrow(model_data),
      n_event = sum(model_data$target_y == 1L),
      n_control = sum(model_data$target_y == 0L),
      stringsAsFactors = FALSE
    )
}

# 使用 10 个 fork 进程并行拟合；错误回传主进程后直接停止
parallel_results <- parallel::mclapply(
  seq_len(nrow(model_tasks)),
  function(task_index) {
    tryCatch(
      list(result = fit_retinal_feature_cox(task_index), error = NULL),
      error = function(error) {
        list(result = NULL, error = conditionMessage(error))
      }
    )
  },
  mc.cores = n_threads,
  mc.preschedule = TRUE
)

failed_tasks <- which(vapply(
  parallel_results,
  function(task_result) !is.null(task_result$error),
  logical(1)
))
if (length(failed_tasks) > 0L) {
  stop(
    "Parallel retinal feature Cox model failed: ",
    paste(
      vapply(
        parallel_results[failed_tasks],
        function(task_result) task_result$error,
        character(1)
      ),
      collapse = " | "
    )
  )
}

all_model_results <- lapply(
  parallel_results,
  function(task_result) task_result$result
)

cox_results <- do.call(rbind, all_model_results)
rownames(cox_results) <- NULL

if (nrow(cox_results) != 290L) {
  stop("Retinal feature Cox result does not contain 290 models")
}

# 计算病种内及全局多重检验校正
cox_results$p_fdr_disease <- ave(
  cox_results$p_value,
  cox_results$disease,
  FUN = function(x) p.adjust(x, method = "BH")
)
cox_results$p_fdr_global <- p.adjust(
  cox_results$p_value,
  method = "BH"
)
cox_results$p_bonferroni_global <- p.adjust(
  cox_results$p_value,
  method = "bonferroni"
)

# 汇总每病种原始 discovery 事件数和逐表型完整病例范围
sample_summary_parts <- lapply(disease_names, function(disease_name) {
  disease_data <- disease_data_list[[disease_name]]
  disease_results <- cox_results[
    cox_results$disease == disease_name,
    ,
    drop = FALSE
  ]

  data.frame(
    disease = disease_name,
    n_discovery_outcome = nrow(disease_data),
    n_event_discovery_outcome = sum(disease_data$target_y == 1L),
    n_control_discovery_outcome = sum(disease_data$target_y == 0L),
    n_models = nrow(disease_results),
    n_model_min = min(disease_results$n),
    n_model_median = median(disease_results$n),
    n_model_max = max(disease_results$n),
    n_event_min = min(disease_results$n_event),
    n_event_median = median(disease_results$n_event),
    n_event_max = max(disease_results$n_event),
    stringsAsFactors = FALSE
  )
})
sample_summary <- do.call(rbind, sample_summary_parts)
rownames(sample_summary) <- NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cox_results, result_path, row.names = FALSE, na = "")
write.csv(sample_summary, sample_summary_path, row.names = FALSE, na = "")

print(sample_summary)
