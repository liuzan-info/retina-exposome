# Step3 independent-exposure selection; original clustering and elimination rules.
library(tidyverse)
library(fs)
library(car)
library(caret)
library(broom)


association_file <- "data/discovery_exposure_associations.csv"
exposure_dic_file <- "data/exposure_dictionary.csv"
oct_dic_file <- "data/retinal_dictionary.csv"
oct_feature_file <- "data/retinal_phenotypes.tsv"
exposure_file <- "data/discovery_exposures.csv"
covariate_file <- "data/gwas_covariates.tsv"
output_root <- "results/independent_exposures"
sig_dir <- file.path(output_root, "01_sig_association")
cluster_vif_dir <- file.path(output_root, "02_cluster_vif")
vif_filtered_dir <- file.path(output_root, "03_vif_threshold3")
pvalue_dir <- file.path(output_root, "04_pvalue5e-2")
summary_dir <- file.path(output_root, "summary")
fs::dir_create(c(sig_dir, cluster_vif_dir, vif_filtered_dir, pvalue_dir, summary_dir))
height_threshold <- 0.3
correlation_cutoff <- 0.90
vif_threshold <- 3
pvalue_threshold <- 5e-2
covars <- c("Age", "Sex", "avMSE")

exposure_dic <- read_csv(exposure_dic_file, col_types = cols(.default = col_character()), show_col_types = FALSE) |>
  dplyr::mutate(semantic_numeric_opposite = as.logical(semantic_numeric_opposite))

association_discovery <- read_csv(association_file, show_col_types = FALSE) |>
  dplyr::semi_join(exposure_dic, by = c("term" = "association_exposure"))
oct_dic <- read_csv(oct_dic_file, show_col_types = FALSE) |>
  dplyr::rename_with(~ str_remove(.x, "^\ufeff"))

sig_p_cutoff <- 0.05 /
  n_distinct(association_discovery$oct_features) /
  n_distinct(association_discovery$term)

sig_associate <- association_discovery |>
  dplyr::transmute(
    EyeFeature = oct_features,
    Exposure = term,
    discovery_estimate = estimate,
    discovery_std.error = std.error,
    discovery_statistic = statistic,
    discovery_p.value = p.value,
    discovery_conf.low = conf.low,
    discovery_conf.high = conf.high,
    discovery_nobs = nobs
  ) |>
  dplyr::mutate(
    discovery_p.value.bonferroni = p.adjust(discovery_p.value, method = "bonferroni"),
    discovery_p.value.fdr = p.adjust(discovery_p.value, method = "BH")
  ) |>
  dplyr::filter(discovery_p.value < sig_p_cutoff) |>
  dplyr::mutate(association_exposure = Exposure, .after = Exposure) |>
  dplyr::left_join(exposure_dic, by = c("Exposure" = "association_exposure")) |>
  dplyr::mutate(
    discovery_estimate_semantic = if_else(
      semantic_numeric_opposite,
      -discovery_estimate,
      discovery_estimate
    )
  ) |>
  dplyr::arrange(EyeFeature, discovery_p.value, Exposure)

write_csv(sig_associate, file.path(sig_dir, "sig_association_discovery.csv"))

selected_exposure_map <- sig_associate |>
  dplyr::transmute(
    expression_exposure,
    association_exposure = Exposure,
    semantic_numeric_opposite
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(association_exposure)

oct_feature_df <- read_table(oct_feature_file, show_col_types = FALSE) |>
  dplyr::rename(eid = IID) |>
  dplyr::select(eid, all_of(unique(sig_associate$EyeFeature)))

covariate_df <- read_table(covariate_file, show_col_types = FALSE) |>
  dplyr::rename(eid = IID) |>
  dplyr::select(eid, all_of(covars)) |>
  dplyr::mutate(Sex = as.factor(Sex))

rename_exposure_map <- setNames(
  selected_exposure_map$expression_exposure,
  selected_exposure_map$association_exposure
)

exposure_df <- read_csv(exposure_file, show_col_types = FALSE) |>
  dplyr::select(eid, all_of(selected_exposure_map$expression_exposure)) |>
  dplyr::rename(!!!rename_exposure_map)

semantic_reverse_exposures <- selected_exposure_map |>
  dplyr::filter(semantic_numeric_opposite) |>
  dplyr::pull(association_exposure)

# 将反向数值编码镜像到统一语义方向后，再用于距离聚类和后续筛选
exposure_df <- exposure_df |>
  dplyr::mutate(
    across(
      all_of(semantic_reverse_exposures),
      ~ min(.x, na.rm = TRUE) + max(.x, na.rm = TRUE) - .x
    )
  )

analysis_df <- oct_feature_df |>
  dplyr::left_join(covariate_df, by = "eid") |>
  dplyr::left_join(exposure_df, by = "eid")

make_lm_formula <- function(response, predictors) {
  as.formula(paste0("`", response, "` ~ ", paste0("`", predictors, "`", collapse = " + ")))
}

complete_n <- function(df, vars) {
  df |>
    dplyr::select(all_of(vars)) |>
    na.omit() |>
    nrow()
}

find_cor_na_exposures <- function(cor_matrix) {
  na_position <- which(is.na(cor_matrix), arr.ind = TRUE)
  na_position <- na_position[na_position[, 1] != na_position[, 2], , drop = FALSE]

  if (nrow(na_position) == 0) {
    return(character())
  }

  unique(c(
    rownames(cor_matrix)[na_position[, 1]],
    colnames(cor_matrix)[na_position[, 2]]
  ))
}

choose_exposure_to_remove <- function(sig_associate, eye_feature, exposures) {
  ranked <- sig_associate |>
    dplyr::filter(EyeFeature == eye_feature, Exposure %in% exposures) |>
    dplyr::arrange(desc(discovery_p.value), discovery_nobs) |>
    dplyr::pull(Exposure)

  if (length(ranked) == 0) exposures[1] else ranked[1]
}

choose_by_discovery_rank <- function(sig_associate, eye_feature, exposures) {
  choose_exposure_to_remove(sig_associate, eye_feature, exposures)
}

calculate_vif_for_cluster <- function(analysis_df, sig_associate, eye_feature, cluster_exposure, correlation_cutoff = 0.90) {
  cluster_nobs <- complete_n(analysis_df, c(eye_feature, cluster_exposure))

  if (length(cluster_exposure) == 1) {
    return(tibble(
      exposure = cluster_exposure,
      vif = 1,
      vif_na_reason = "",
      vif_nobs = cluster_nobs
    ))
  }

  tryCatch({
    kept_for_correlation <- cluster_exposure
    cor_na_removed_exposure <- character()

    while (length(kept_for_correlation) > 1) {
      cluster_values <- analysis_df |>
        dplyr::select(all_of(kept_for_correlation))
      cor_matrix <- cor(cluster_values, use = "pairwise.complete.obs")

      if (!any(is.na(cor_matrix))) {
        break
      }

      na_exposures <- find_cor_na_exposures(cor_matrix)
      if (length(na_exposures) == 0) {
        na_exposures <- kept_for_correlation
      }
      remove_one <- choose_exposure_to_remove(sig_associate, eye_feature, na_exposures)
      kept_for_correlation <- setdiff(kept_for_correlation, remove_one)
      cor_na_removed_exposure <- c(cor_na_removed_exposure, remove_one)
    }

    cor_na_removed_df <- tibble(exposure = cor_na_removed_exposure) |>
      dplyr::mutate(
        vif = NA_real_,
        vif_na_reason = "Removed due to NA in correlation matrix",
        vif_nobs = cluster_nobs
      )

    if (length(kept_for_correlation) == 1) {
      kept_df <- tibble(
        exposure = kept_for_correlation,
        vif = 1,
        vif_na_reason = "",
        vif_nobs = complete_n(analysis_df, c(eye_feature, kept_for_correlation))
      )
      return(dplyr::bind_rows(kept_df, cor_na_removed_df))
    }

    cluster_values <- analysis_df |>
      dplyr::select(all_of(kept_for_correlation))
    cor_matrix <- cor(cluster_values, use = "pairwise.complete.obs")

    if (any(is.na(cor_matrix))) {
      stop("Correlation matrix still contains NA after removing NA-triggering exposures")
    }

    remove_idx <- caret::findCorrelation(cor_matrix, cutoff = correlation_cutoff)
    removed_exposure <- if (length(remove_idx) == 0) character() else colnames(cor_matrix)[remove_idx]
    kept_exposure <- setdiff(kept_for_correlation, removed_exposure)

    removed_df <- tibble(
      exposure = removed_exposure,
      vif = NA_real_,
      vif_na_reason = "Removed due to collinearity (correlation > 0.90)",
      vif_nobs = cluster_nobs
    ) |>
      dplyr::bind_rows(cor_na_removed_df)

    if (length(kept_exposure) == 1) {
      kept_df <- tibble(
        exposure = kept_exposure,
        vif = 1,
        vif_na_reason = "",
        vif_nobs = complete_n(analysis_df, c(eye_feature, kept_exposure))
      )
      return(dplyr::bind_rows(kept_df, removed_df))
    }

    model_data <- analysis_df |>
      dplyr::select(all_of(c(eye_feature, kept_exposure)))
    model <- lm(make_lm_formula(eye_feature, kept_exposure), data = model_data)
    vif_result <- car::vif(model)
    names(vif_result) <- str_remove_all(names(vif_result), "`")

    kept_df <- tibble(
      exposure = names(vif_result),
      vif = as.numeric(vif_result),
      vif_na_reason = "",
      vif_nobs = nobs(model)
    )

    dplyr::bind_rows(kept_df, removed_df)
  }, error = function(e) {
    tibble(
      exposure = cluster_exposure,
      vif = NA_real_,
      vif_na_reason = paste0("VIF calculation error: ", e$message),
      vif_nobs = NA_integer_
    )
  })
}

recalculate_vif <- function(analysis_df, eye_feature, exposures) {
  if (length(exposures) == 1) {
    return(tibble(exposure = exposures, vif = 1))
  }

  model_data <- analysis_df |>
    dplyr::select(all_of(c(eye_feature, exposures)))
  model <- lm(make_lm_formula(eye_feature, exposures), data = model_data)
  vif_result <- car::vif(model)
  names(vif_result) <- str_remove_all(names(vif_result), "`")

  tibble(
    exposure = names(vif_result),
    vif = as.numeric(vif_result)
  )
}

calculate_exposure_pvalues <- function(analysis_df, eye_feature, exposures, covars) {
  tryCatch({
    model_data <- analysis_df |>
      dplyr::select(all_of(c(eye_feature, covars, exposures)))
    model <- lm(make_lm_formula(eye_feature, c(covars, exposures)), data = model_data)
    model_result <- broom::tidy(model)

    model_result |>
      dplyr::filter(term %in% exposures) |>
      dplyr::transmute(
        exposure = term,
        pvalue = p.value,
        pvalue_na_reason = "",
        pvalue_nobs = nobs(model)
      )
  }, error = function(e) {
    tibble(
      exposure = exposures,
      pvalue = NA_real_,
      pvalue_na_reason = paste0("P-value model error: ", e$message),
      pvalue_nobs = NA_integer_
    )
  })
}

find_key_exposure <- function(analysis_df, sig_associate, eye_feature, output_dir, height_threshold = 0.3, correlation_cutoff = 0.90) {
  eye_output_dir <- file.path(output_dir, eye_feature)
  dir_create(eye_output_dir)

  eye_exposures <- sig_associate |>
    dplyr::filter(EyeFeature == eye_feature) |>
    dplyr::pull(Exposure) |>
    unique()

  eye_exposure_data <- analysis_df |>
    dplyr::select(all_of(eye_exposures))

  dist_removed_exposures <- character()

  if (length(eye_exposures) == 1) {
    cluster_df <- tibble(exposure = eye_exposures, cluster_id = 1L)
  } else {
    while (ncol(eye_exposure_data) > 1) {
      dist_matrix <- dist(t(eye_exposure_data), method = "euclidean")
      if (!any(is.na(dist_matrix))) break

      na_position <- which(is.na(as.matrix(dist_matrix)), arr.ind = TRUE)
      na_position <- unique(t(apply(na_position, 1, sort)))
      na_exposure <- unique(c(
        rownames(as.matrix(dist_matrix))[na_position[, 1]],
        colnames(as.matrix(dist_matrix))[na_position[, 2]]
      ))

      remove_one <- choose_by_discovery_rank(sig_associate, eye_feature, na_exposure)
      eye_exposure_data <- eye_exposure_data |>
        dplyr::select(-all_of(remove_one))
      dist_removed_exposures <- c(dist_removed_exposures, remove_one)
    }

    if (ncol(eye_exposure_data) == 1) {
      cluster_df <- tibble(exposure = colnames(eye_exposure_data), cluster_id = 1L)
    } else {
      dist_matrix <- dist(t(eye_exposure_data), method = "euclidean")
      hc <- hclust(dist_matrix, method = "ward.D2")
      n_exposure <- ncol(eye_exposure_data)
      cut_index <- n_exposure - ceiling(n_exposure * height_threshold)
      cut_index <- max(1, min(length(hc$height), cut_index))
      cut_height <- hc$height[cut_index]
      exposure_clusters <- cutree(hc, h = cut_height)

      cluster_df <- tibble(
        exposure = names(exposure_clusters),
        cluster_id = as.integer(exposure_clusters)
      )
    }
  }

  vif_df <- map_dfr(sort(unique(cluster_df$cluster_id)), function(cluster_i) {
    cluster_exposure <- cluster_df |>
      dplyr::filter(cluster_id == cluster_i) |>
      dplyr::pull(exposure)

    calculate_vif_for_cluster(
      analysis_df = analysis_df,
      sig_associate = sig_associate,
      eye_feature = eye_feature,
      cluster_exposure = cluster_exposure,
      correlation_cutoff = correlation_cutoff
    )
  })

  cluster_df <- cluster_df |>
    dplyr::left_join(vif_df, by = "exposure") |>
    dplyr::arrange(cluster_id, exposure)

  if (length(dist_removed_exposures) > 0) {
    cluster_df <- dplyr::bind_rows(
      cluster_df,
      tibble(
        exposure = dist_removed_exposures,
        cluster_id = NA_integer_,
        vif = NA_real_,
        vif_na_reason = "Removed due to distance calculation issue (NA in dist_matrix)",
        vif_nobs = NA_integer_
      )
    )
  }

  write_csv(cluster_df, file.path(eye_output_dir, "exposure_clusterId_vif.csv"))
}

for (eye_feature in unique(sig_associate$EyeFeature)) {
  message(paste0("Cluster and VIF: ", eye_feature))
  find_key_exposure(
    analysis_df = analysis_df,
    sig_associate = sig_associate,
    eye_feature = eye_feature,
    output_dir = cluster_vif_dir,
    height_threshold = height_threshold,
    correlation_cutoff = correlation_cutoff
  )
}

remove_high_vif <- function(analysis_df, sig_associate, eye_feature, input_dir, output_dir, vif_threshold = 3) {
  eye_output_dir <- file.path(output_dir, eye_feature)
  dir_create(eye_output_dir)

  exposure_cluster_vif <- read_csv(
    file.path(input_dir, eye_feature, "exposure_clusterId_vif.csv"),
    col_types = cols(exposure = col_character(), .default = col_guess())
  ) |>
    dplyr::filter(!is.na(cluster_id), !is.na(vif))

  filtered_clusters <- map_dfr(sort(unique(exposure_cluster_vif$cluster_id)), function(cluster_i) {
    current_cluster <- exposure_cluster_vif |>
      dplyr::filter(cluster_id == cluster_i)

    while (any(current_cluster$vif > vif_threshold, na.rm = TRUE)) {
      max_vif <- max(current_cluster$vif, na.rm = TRUE)
      candidate_exposures <- current_cluster |>
        dplyr::filter(vif == max_vif) |>
        dplyr::pull(exposure)

      remove_one <- choose_by_discovery_rank(sig_associate, eye_feature, candidate_exposures)
      current_cluster <- current_cluster |>
        dplyr::filter(exposure != remove_one)

      if (nrow(current_cluster) == 1) {
        current_cluster$vif <- 1
      } else {
        updated_vif <- tryCatch(
          recalculate_vif(analysis_df, eye_feature, current_cluster$exposure),
          error = function(e) tibble(exposure = current_cluster$exposure, vif = NA_real_)
        )

        current_cluster <- current_cluster |>
          dplyr::select(-vif) |>
          dplyr::left_join(updated_vif, by = "exposure")
      }
    }

    current_cluster
  })

  write_csv(filtered_clusters, file.path(eye_output_dir, "exposure_clusterId_vif_filtered.csv"))
}

for (eye_feature in unique(sig_associate$EyeFeature)) {
  message(paste0("Remove high VIF: ", eye_feature))
  remove_high_vif(
    analysis_df = analysis_df,
    sig_associate = sig_associate,
    eye_feature = eye_feature,
    input_dir = cluster_vif_dir,
    output_dir = vif_filtered_dir,
    vif_threshold = vif_threshold
  )
}

remove_high_pvalue <- function(analysis_df, eye_feature, input_dir, output_dir, pvalue_threshold, covars = c("Age", "Sex", "avMSE")) {
  eye_output_dir <- file.path(output_dir, eye_feature)
  dir_create(eye_output_dir)

  exposure_cluster_vif <- read_csv(
    file.path(input_dir, eye_feature, "exposure_clusterId_vif_filtered.csv"),
    col_types = cols(exposure = col_character(), .default = col_guess())
  )

  filtered_clusters <- map_dfr(sort(unique(exposure_cluster_vif$cluster_id)), function(cluster_i) {
    current_cluster <- exposure_cluster_vif |>
      dplyr::filter(cluster_id == cluster_i)

    current_cluster$pvalue <- NA_real_
    current_cluster$pvalue_na_reason <- ""
    current_cluster$pvalue_nobs <- NA_integer_

    repeat {
      if (nrow(current_cluster) <= 1) break

      pvalue_df <- calculate_exposure_pvalues(
        analysis_df = analysis_df,
        eye_feature = eye_feature,
        exposures = current_cluster$exposure,
        covars = covars
      )

      current_cluster <- current_cluster |>
        dplyr::select(-pvalue, -pvalue_na_reason, -pvalue_nobs) |>
        dplyr::left_join(pvalue_df, by = "exposure")

      if (any(is.na(current_cluster$pvalue))) {
        current_cluster <- current_cluster |>
          dplyr::filter(!is.na(pvalue))
        next
      }

      high_pvalue <- current_cluster |>
        dplyr::filter(pvalue > pvalue_threshold)

      if (nrow(high_pvalue) == 0) break

      remove_one <- high_pvalue |>
        dplyr::arrange(desc(pvalue)) |>
        dplyr::slice_head(n = 1) |>
        dplyr::pull(exposure)

      current_cluster <- current_cluster |>
        dplyr::filter(exposure != remove_one)
    }

    current_cluster
  })

  write_csv(filtered_clusters, file.path(eye_output_dir, "exposure_clusterId_vif_pvalue.csv"))
}

for (eye_feature in unique(sig_associate$EyeFeature)) {
  message(paste0("Remove high p-value: ", eye_feature))
  remove_high_pvalue(
    analysis_df = analysis_df,
    eye_feature = eye_feature,
    input_dir = vif_filtered_dir,
    output_dir = pvalue_dir,
    pvalue_threshold = pvalue_threshold,
    covars = covars
  )
}

final_result_files <- dir_ls(
  pvalue_dir,
  recurse = TRUE,
  regexp = "exposure_clusterId_vif_pvalue\\.csv$"
)

key_exposures <- map_dfr(final_result_files, function(file) {
  read_csv(file, col_types = cols(exposure = col_character(), .default = col_guess())) |>
    dplyr::mutate(EyeFeature = path_file(path_dir(file)), .before = 1)
})

write_csv(key_exposures, file.path(summary_dir, "key_exposures_pvalue5e-2.csv"))

key_exposures_annotated <- key_exposures |>
  dplyr::left_join(oct_dic, by = "EyeFeature") |>
  dplyr::left_join(
    exposure_dic |>
      dplyr::select(
        expression_exposure,
        association_exposure,
        exposure_modifiability,
        semantic_numeric_opposite,
        exposure_value_meaning_simple,
        exposure_plot_label,
        exposure_value_meaning,
        Category,
        FieldID,
        Field,
        Type,
        Sexed
      ),
    by = c("exposure" = "association_exposure")
  ) |>
  dplyr::left_join(
    sig_associate |>
      dplyr::select(
        EyeFeature,
        Exposure,
        discovery_estimate,
        discovery_estimate_semantic,
        discovery_p.value,
        discovery_nobs,
        discovery_p.value.fdr,
        discovery_p.value.bonferroni
      ),
    by = c("EyeFeature", "exposure" = "Exposure")
  )

write_csv(key_exposures_annotated, file.path(summary_dir, "key_exposures_pvalue5e-2_annotated.csv"))
