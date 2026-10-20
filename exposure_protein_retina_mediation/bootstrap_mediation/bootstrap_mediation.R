# Bootstrap protein mediation on the fixed discovery+replication candidate paths.
library(data.table)
library(mediation)
library(foreach)
library(doParallel)
PATH_TRIPLES <- "data/mediation_candidates.csv"
PATH_PROTEIN <- "data/pooled_mediation_proteins.csv"
PATH_COV_DISCOVERY <- "data/discovery_protein_covariates.csv"
PATH_COV_REP1 <- "data/rep1_protein_covariates.csv"
PATH_EXPOSURE <- "data/pooled_exposures.csv"
PATH_EXPOSURE_DICTIONARY <- "data/exposure_dictionary.csv"
PATH_OCT <- "data/retinal_phenotypes.tsv"
PATH_COHORT <- "data/pooled_keep.tsv"
PATH_OUTPUT_DIR <- "results/protein_mediation"
PATH_RESULT <- "results/protein_mediation/mediation_results_boot2000.csv"
N_THREADS <- 50L
SIMS <- 2000L
MIN_N <- 50L
CHUNK_SIZE <- 100L
COVARIATES <- c("Age", "Sex", "avMSE", "eGFR")

set.seed(20260721)
dir.create(PATH_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 读取候选三元组并构建 discovery+rep1 合并分析数据
triples <- fread(PATH_TRIPLES)

if (!all(c("exposure", "protein", "EyeFeature") %in% names(triples))) {
  stop("三元组文件缺少 exposure、protein 或 EyeFeature 列。")
}
if (nrow(triples) != 4528L) {
  stop(sprintf("三元组数应为 4528，当前为 %d。", nrow(triples)))
}
if (anyDuplicated(triples, by = c("exposure", "protein", "EyeFeature"))) {
  stop("三元组三键存在重复。")
}

triples[, task_id := .I]
exposure_vars <- unique(triples$exposure)
protein_vars <- unique(triples$protein)
oct_vars <- unique(triples$EyeFeature)

if (length(exposure_vars) != 67L || length(protein_vars) != 177L || length(oct_vars) != 32L) {
  stop("候选三元组的暴露、蛋白或 OCT 特征数量与预期不一致。")
}

exposure_dictionary <- fread(
  PATH_EXPOSURE_DICTIONARY,
  select = c("association_exposure", "expression_exposure")
)
if (anyDuplicated(exposure_dictionary, by = "association_exposure")) {
  stop("最终暴露词典的 association_exposure 存在重复。")
}

exposure_map <- merge(
  data.table(exposure = exposure_vars),
  exposure_dictionary,
  by.x = "exposure",
  by.y = "association_exposure",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(exposure_map$expression_exposure)) {
  stop("候选关联暴露无法映射到 expression_exposure。")
}

expression_exposure_vars <- exposure_map$expression_exposure
exposure_header <- names(fread(PATH_EXPOSURE, nrows = 0L))
missing_expression_header <- setdiff(expression_exposure_vars, exposure_header)
if (length(missing_expression_header) > 0L) {
  stop("暴露矩阵表头缺少映射后的 expression_exposure。")
}

exposure_map[, X := sprintf("X_%03d", seq_len(.N))]

cohort <- fread(PATH_COHORT, select = "IID")
setnames(cohort, "IID", "eid")
cohort[, eid := as.character(eid)]

exposure_data <- fread(PATH_EXPOSURE, select = c("eid", expression_exposure_vars))
protein_data <- fread(PATH_PROTEIN, select = c("eid", protein_vars))
oct_data <- fread(PATH_OCT, select = c("IID", oct_vars))
cov_discovery <- fread(PATH_COV_DISCOVERY, select = c("eid", COVARIATES))
cov_rep1 <- fread(PATH_COV_REP1, select = c("eid", COVARIATES))

setnames(oct_data, "IID", "eid")
covariate_data <- rbindlist(list(cov_discovery, cov_rep1), use.names = TRUE)

for (data_name in c("exposure_data", "protein_data", "oct_data", "covariate_data")) {
  get(data_name)[, eid := as.character(eid)]
}

for (data_name in c("exposure_data", "protein_data", "oct_data", "covariate_data")) {
  assign(data_name, merge(cohort, get(data_name), by = "eid"))
}

missing_expression_exposure <- setdiff(expression_exposure_vars, names(exposure_data))
missing_protein <- setdiff(protein_vars, names(protein_data))
missing_oct <- setdiff(oct_vars, names(oct_data))
missing_covariates <- setdiff(COVARIATES, names(covariate_data))
if (length(missing_expression_exposure) > 0L || length(missing_protein) > 0L ||
    length(missing_oct) > 0L || length(missing_covariates) > 0L) {
  stop("输入数据缺少候选变量或协变量。")
}

protein_map <- data.table(protein = protein_vars, M = sprintf("M_%03d", seq_along(protein_vars)))
oct_map <- data.table(EyeFeature = oct_vars, Y = sprintf("Y_%03d", seq_along(oct_vars)))

setnames(exposure_data, exposure_map$expression_exposure, exposure_map$X)
setnames(protein_data, protein_map$protein, protein_map$M)
setnames(oct_data, oct_map$EyeFeature, oct_map$Y)

full_data <- merge(exposure_data, protein_data, by = "eid", all.x = TRUE)
full_data <- merge(full_data, oct_data, by = "eid", all.x = TRUE)
full_data <- merge(full_data, covariate_data, by = "eid", all.x = TRUE)
full_data[, Sex := as.factor(Sex)]

tasks <- merge(triples, exposure_map, by = "exposure", all.x = TRUE, sort = FALSE)
tasks <- merge(tasks, protein_map, by = "protein", all.x = TRUE, sort = FALSE)
tasks <- merge(tasks, oct_map, by = "EyeFeature", all.x = TRUE, sort = FALSE)
setorder(tasks, task_id)

if (anyNA(tasks$X) || anyNA(tasks$M) || anyNA(tasks$Y)) {
  stop("候选三元组无法映射到内部变量名。")
}
if (!all(unique(c(tasks$X, tasks$M, tasks$Y, COVARIATES)) %in% names(full_data))) {
  stop("合并分析数据缺少建模变量。")
}

cat(sprintf("候选三元组：%d；暴露：%d；蛋白：%d；OCT 特征：%d。\n",
            nrow(tasks), length(exposure_vars), length(protein_vars), length(oct_vars)))
cat(sprintf("合并后的 discovery+rep1 分析数据：%d 人。\n", nrow(full_data)))

n_tasks <- nrow(tasks)
n_chunks <- ceiling(n_tasks / CHUNK_SIZE)
cl <- parallel::makeCluster(N_THREADS)
parallel::clusterCall(cl, function() { library(data.table); library(mediation); NULL })
doParallel::registerDoParallel(cl)
parallel::clusterSetRNGStream(cl, iseed = 20260721)
chunks <- vector("list", n_chunks)
analysis_start <- Sys.time()
for (chunk_idx in seq_len(n_chunks)) {
  start_idx <- (chunk_idx - 1L) * CHUNK_SIZE + 1L
  end_idx <- min(chunk_idx * CHUNK_SIZE, n_tasks)
    tasks_chunk <- tasks[start_idx:end_idx]
    result_chunk <- foreach(
      i = seq_len(nrow(tasks_chunk)),
      .combine = "rbind"
    ) %dopar% {
      task <- tasks_chunk[i]
      X_var <- task$X
      M_var <- task$M
      Y_var <- task$Y
      required_vars <- c(X_var, M_var, Y_var, COVARIATES)

      analysis_data <- full_data[, ..required_vars]
      analysis_data <- analysis_data[stats::complete.cases(analysis_data)]
      analysis_data <- as.data.frame(analysis_data)

      if (nrow(analysis_data) < MIN_N) {
        stop(sprintf("任务 %d 的完整案例数少于 %d。", task$task_id, MIN_N))
      }
      if (length(unique(analysis_data[[X_var]])) < 2L ||
          length(unique(analysis_data[[M_var]])) < 2L ||
          length(unique(analysis_data[[Y_var]])) < 2L) {
        stop(sprintf("任务 %d 的暴露、蛋白或 OCT 特征无变异。", task$task_id))
      }

      mediator_formula <- stats::reformulate(c(X_var, COVARIATES), response = M_var)
      outcome_formula <- stats::reformulate(c(X_var, M_var, COVARIATES), response = Y_var)
      mediator_fit <- stats::lm(mediator_formula, data = analysis_data)
      outcome_fit <- stats::lm(outcome_formula, data = analysis_data)
      mediator_fit$call$formula <- mediator_formula
      outcome_fit$call$formula <- outcome_formula

      mediation_fit <- mediation::mediate(
        model.m = mediator_fit,
        model.y = outcome_fit,
        treat = X_var,
        mediator = M_var,
        boot = TRUE,
        sims = SIMS
      )

      task_output <- as.data.frame(task[, setdiff(names(task), c("X", "M", "Y")), with = FALSE])
      cbind(
        task_output,
        data.frame(
          n_obs = nrow(analysis_data),
          acme_est = mediation_fit$d.avg,
          acme_p = mediation_fit$d.avg.p,
          acme_ci_low = mediation_fit$d.avg.ci[1],
          acme_ci_high = mediation_fit$d.avg.ci[2],
          ade_est = mediation_fit$z.avg,
          ade_p = mediation_fit$z.avg.p,
          ade_ci_low = mediation_fit$z.avg.ci[1],
          ade_ci_high = mediation_fit$z.avg.ci[2],
          total_est = mediation_fit$tau.coef,
          total_p = mediation_fit$tau.p,
          total_ci_low = mediation_fit$tau.ci[1],
          total_ci_high = mediation_fit$tau.ci[2],
          prop_med_est = mediation_fit$n.avg,
          prop_med_ci_low = mediation_fit$n.avg.ci[1],
          prop_med_ci_high = mediation_fit$n.avg.ci[2],
          sims = SIMS
        )
      )
    }


  chunks[[chunk_idx]] <- result_chunk
  elapsed <- as.numeric(difftime(Sys.time(), analysis_start, units = "secs"))
  message(sprintf("Mediation: %d/%d (%.1f%%); ETA %.1f s", end_idx, n_tasks,
                  100 * end_idx / n_tasks, elapsed / end_idx * (n_tasks - end_idx)))
}
parallel::stopCluster(cl)
results <- data.table::rbindlist(chunks, use.names = TRUE)
data.table::setorder(results, task_id)
data.table::fwrite(results, PATH_RESULT)
