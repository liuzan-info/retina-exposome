# Main retinal-PRS Cox model; full-scoring-population PRS scaling is done before exclusions.
prs_path <- "data/white_non_discovery_prs.tsv"
covariate_path <- "data/disease_covariates.csv"
pca_path <- "data/disease_pca.txt"
dictionary_path <- "data/retinal_dictionary.csv"
disease_dir <- "data/incident_diseases"
output_dir <- "results/retinal_prs_disease"
result_path <- file.path(output_dir, "cox_prs_eye_disease_results.csv")
sample_summary_path <- file.path(output_dir, "sample_summary.csv")
disease_names <- c("AMD", "RD", "Cataract", "POAG", "DR")
n_workers <- 5L

library(data.table)
library(survival)

# 读取 58 个表型名称及注释
feature_dictionary <- data.table::fread(
  dictionary_path,
  select = c("EyeFeature", "feature_name", "plot_label_med", "group")
)
trait_names <- feature_dictionary$EyeFeature
prs_columns <- paste0(trait_names, "_AVG")
if (length(trait_names) != 58L) {
  stop("Retinal phenotype dictionary must contain 58 traits")
}

# 读取白人非 Discovery 人群的 PRS，并按各表型 PRS 队列分布转换为 z 分数
prs_data <- data.table::fread(prs_path, select = c("IID", prs_columns))
prs_data[, IID := as.character(IID)]
if (anyDuplicated(prs_data$IID)) {
  stop("PRS matrix contains duplicated IID")
}
prs_data[, (prs_columns) := lapply(
  .SD,
  function(x) as.numeric(scale(as.numeric(x)))
), .SDcols = prs_columns]

# 按 IID/eid 匹配协变量和前 10 个遗传主成分
covariates <- data.table::fread(
  covariate_path,
  select = c("eid", "Age", "Sex", "avMSE", "gbatch")
)
covariates[, eid := as.character(eid)]
if (anyDuplicated(covariates$eid)) {
  stop("Covariate file contains duplicated eid")
}
pca_columns <- paste0("PC", 1:10)
pca_data <- data.table::as.data.table(utils::read.table(
  pca_path,
  header = FALSE,
  sep = "",
  fill = TRUE,
  col.names = c("eid", paste0("PC", 1:40)),
  colClasses = c("character", rep("numeric", 10), rep("NULL", 30))
))
if (anyDuplicated(pca_data$eid)) {
  stop("PCA file contains duplicated eid")
}
aligned_covariates <- covariates[
  match(prs_data$IID, covariates$eid),
  .(Age, Sex, avMSE, gbatch)
]
aligned_pca <- pca_data[match(prs_data$IID, pca_data$eid), ..pca_columns]
prs_data[, c("Age", "Sex", "avMSE", "gbatch") := aligned_covariates]
prs_data[, (pca_columns) := aligned_pca]
# 正式模型共同排除年龄、性别、分型批次或任一 PC 缺失者
complete_pca <- stats::complete.cases(prs_data[, ..pca_columns])
prs_data <- prs_data[
  !is.na(Age) & !is.na(Sex) & !is.na(gbatch) & complete_pca
]

# 在确定的分析队列内按性别计算中位数，一次填补供五种眼病共用
prs_data[, avMSE_was_missing := is.na(avMSE)]
prs_data[, avMSE_imputed := {
  sex_median <- stats::median(avMSE, na.rm = TRUE)
  replace(avMSE, is.na(avMSE), sex_median)
}, by = Sex]
if (anyNA(prs_data$avMSE_imputed)) {
  stop("Sex-specific avMSE median imputation left missing values")
}

# 提取每个 Cox 模型的 PRS 效应与样本信息
fit_prs_cox <- function(model_data, formula, model_name, disease_name, trait_index) {
  fit <- survival::coxph(
    formula,
    data = model_data,
    ties = "efron",
    na.action = na.fail,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  estimate <- summary(fit)$coefficients["PRS", ]
  beta <- unname(estimate["coef"])
  standard_error <- unname(estimate["se(coef)"])
  feature <- feature_dictionary[trait_index]

  data.table::data.table(
    disease = disease_name,
    retinal_feature = trait_names[[trait_index]],
    prs_column = prs_columns[[trait_index]],
    feature_name = feature$feature_name,
    plot_label = feature$plot_label_med,
    feature_group = feature$group,
    model = model_name,
    effect_unit = "per 1 SD of PRS in the scored cohort",
    beta = beta,
    standard_error = standard_error,
    hazard_ratio = exp(beta),
    confidence_interval_low = exp(beta - stats::qnorm(0.975) * standard_error),
    confidence_interval_high = exp(beta + stats::qnorm(0.975) * standard_error),
    z_value = unname(estimate["z"]),
    p_value = unname(estimate["Pr(>|z|)"]),
    n = nrow(model_data),
    n_event = sum(model_data$event == 1L),
    n_control = sum(model_data$event == 0L),
    n_avMSE_imputed = if (model_name == "age_sex_avMSE_imputed_gbatch_pc10") {
      sum(model_data$avMSE_was_missing)
    } else {
      0L
    }
  )
}

# 每个进程读取一种 incident 眼病，依次拟合 58 个 PRS 的正式模型
fit_disease <- function(disease_name) {
  disease_data <- data.table::fread(
    file.path(disease_dir, paste0(disease_name, ".csv")),
    select = c("eid", "BL2Target_yrs", "target_y")
  )
  disease_data[, eid := as.character(eid)]
  if (anyDuplicated(disease_data$eid)) {
    stop(disease_name, " incident disease file contains duplicated eid")
  }

  disease_index <- match(prs_data$IID, disease_data$eid)
  prs_index <- which(!is.na(disease_index))
  outcome_index <- disease_index[prs_index]

  model_data <- data.table::data.table(
    time = disease_data$BL2Target_yrs[outcome_index],
    event = disease_data$target_y[outcome_index],
    Age = prs_data$Age[prs_index],
    Sex = factor(prs_data$Sex[prs_index]),
    gbatch = factor(prs_data$gbatch[prs_index]),
    avMSE = prs_data$avMSE[prs_index],
    avMSE_imputed = prs_data$avMSE_imputed[prs_index],
    avMSE_was_missing = prs_data$avMSE_was_missing[prs_index]
  )
  model_data[, (pca_columns) := prs_data[prs_index, ..pca_columns]]

  disease_results <- vector("list", length(trait_names))
  result_index <- 1L
  for (trait_index in seq_along(trait_names)) {
    model_data[, PRS := prs_data[[prs_columns[[trait_index]]]][prs_index]]
    basic_complete <- stats::complete.cases(
      model_data[, .(time, event, Age, Sex, PRS)]
    )
    avMSE_complete <- basic_complete & !is.na(model_data$avMSE)
    basic_data <- model_data[basic_complete]
    complete_avMSE_data <- model_data[avMSE_complete]

    disease_results[[result_index]] <- fit_prs_cox(
      basic_data,
      survival::Surv(time, event) ~ PRS + Age + Sex + avMSE_imputed + gbatch +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
      "age_sex_avMSE_imputed_gbatch_pc10",
      disease_name,
      trait_index
    )
    result_index <- result_index + 1L
  }

  data.table::rbindlist(disease_results)
}

# fork 进程共享已读取的 PRS 数据，每种疾病由一个进程处理
parallel_results <- parallel::mclapply(
  disease_names,
  fit_disease,
  mc.cores = n_workers,
  mc.preschedule = TRUE
)
cox_results <- data.table::rbindlist(parallel_results)
if (nrow(cox_results) != 58L * length(disease_names)) {
  stop("Cox result count is not 58 traits x 5 diseases")
}

# 正式模型分别按病种内 58 项及全局 290 项进行校正
cox_results[, p_fdr_disease := stats::p.adjust(p_value, method = "BH"),
            by = .(model, disease)]
cox_results[, p_fdr_global := stats::p.adjust(p_value, method = "BH"),
            by = model]
cox_results[, p_bonferroni_global := stats::p.adjust(
  p_value, method = "bonferroni"
), by = model]

# 汇总每种疾病和模型的实际人数、事件数及填补人数范围
sample_summary <- cox_results[, .(
  n_models = .N,
  n_min = min(n),
  n_median = stats::median(n),
  n_max = max(n),
  n_event_min = min(n_event),
  n_event_median = stats::median(n_event),
  n_event_max = max(n_event),
  n_avMSE_imputed_min = min(n_avMSE_imputed),
  n_avMSE_imputed_max = max(n_avMSE_imputed)
), by = .(disease, model)]

# 保存两张 CSV；不输出 HTML 以外的图件
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(cox_results, result_path)
data.table::fwrite(sample_summary, sample_summary_path)
