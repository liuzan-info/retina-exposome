# GO BP ORA followed by semantic reduction; inputs are already mapped gene lists.
foreground <- readr::read_csv("data/mediator_foreground_genes.csv", show_col_types = FALSE)$gene
background <- readr::read_csv("data/measured_protein_go_background.csv", show_col_types = FALSE)$gene
msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
stopifnot(identical(unique(msig$db_version), "2025.1.Hs"))
term_gene <- msig |> dplyr::distinct(gs_name, gene_symbol) |>
  dplyr::filter(gene_symbol %in% background)
term_size <- term_gene |> dplyr::count(gs_name, name = "pathway_size") |>
  dplyr::filter(pathway_size >= 10L, pathway_size <= 500L)
term_gene <- term_gene |> dplyr::semi_join(term_size, by = "gs_name")
pathways <- split(term_gene$gene_symbol, term_gene$gs_name)
ora <- dplyr::bind_rows(lapply(names(pathways), function(term) {
  members <- pathways[[term]]
  hits <- intersect(foreground, members)
  k <- length(hits)
  tibble::tibble(pathway_id = term, overlap = k, foreground_size = length(foreground),
    pathway_size = length(members), background_size = length(background),
    fold_enrichment = (k / length(foreground)) / (length(members) / length(background)),
    p_value = stats::phyper(k - 1L, length(members), length(background) - length(members),
                           length(foreground), lower.tail = FALSE),
    overlap_genes = paste(sort(hits), collapse = ";"))
})) |>
  dplyr::mutate(fdr_bh = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::left_join(msig |> dplyr::distinct(gs_name, gs_exact_source, gs_description),
                   by = c("pathway_id" = "gs_name")) |>
  dplyr::arrange(fdr_bh, p_value, pathway_id)
go <- ora |> dplyr::filter(fdr_bh < 0.05)
scores <- stats::setNames(-log10(go$fdr_bh), go$gs_exact_source)
sim_matrix <- rrvgo::calculateSimMatrix(names(scores), orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
reduced <- rrvgo::reduceSimMatrix(sim_matrix, scores, threshold = 0.8, orgdb = "org.Hs.eg.db")
representative_ids <- reduced |> dplyr::filter(parentTerm == term) |> dplyr::pull(go)
mapping <- go |>
  dplyr::left_join(reduced |> dplyr::transmute(go_id = go, cluster, go_name = term,
    representative_go_id = parent, representative_go_name = parentTerm, score),
    by = c("gs_exact_source" = "go_id")) |>
  dplyr::mutate(is_representative = gs_exact_source %in% representative_ids)
representatives <- mapping |> dplyr::filter(is_representative) |>
  dplyr::arrange(fdr_bh, p_value, gs_exact_source)
dir.create("results/go_enrichment", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(ora, "results/go_enrichment/go_bp_complete.csv")
readr::write_csv(mapping, "results/go_enrichment/go_bp_term_mapping.csv")
readr::write_csv(representatives, "results/go_enrichment/go_bp_representative_terms.csv")
