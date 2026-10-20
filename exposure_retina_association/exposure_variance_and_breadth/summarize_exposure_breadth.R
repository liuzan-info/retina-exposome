# Breadth across 19 representative traits and median Delta R-squared among significant pairs.
library(data.table)
pairs <- data.table::fread("results/exposure_variance/pair_results.tsv")
footprint <- pairs[, .(breadth = sum(significant),
  associated_phenotype_ids = paste(EyeFeature[significant], collapse = "; ")),
  by = exposure]
magnitude <- pairs[significant == TRUE,
  .(median_delta_r2_pct = median(delta_r2_pct)), by = exposure]
footprint <- merge(footprint, magnitude, by = "exposure", all.x = TRUE, sort = FALSE)
data.table::fwrite(footprint, "results/exposure_variance/exposure_footprint.tsv", sep = "\t", na = "NA")
