# =========================================================
# Task 4: GO BP enrichment for species-specific and conserved OCRs
# Using rGREAT with species-matched OCR universes as background
# =========================================================

# -----------------------------
# 0. Install packages if needed
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("rtracklayer", "rGREAT", "GenomicRanges", "GenomeInfoDb")
cran_pkgs <- c("ggplot2", "dplyr", "stringr")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(rGREAT)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------
# 1. User settings
# -----------------------------
mouse_tss_source <- "mm10"
human_tss_source <- "hg38"
ontology_to_run  <- "GO:BP"

input_dir <- Sys.getenv("TASK4_INPUT_DIR", unset = "results/mapping")
outdir    <- Sys.getenv("TASK4_OUTDIR", unset = "results/task_4_go_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Input files
# -----------------------------
mouse_specific_file <- file.path(input_dir, "mouse_specific.bed")
human_specific_file <- file.path(input_dir, "human_specific.bed")
conserved_hm_file   <- file.path(input_dir, "conserved_human_in_mouse.bed")

mouse_bg_file <- file.path(input_dir, "mouse_pancreas_ocr.processed.bed")
human_bg_file <- file.path(input_dir, "human_pancreas_ocr.processed.bed")


# -----------------------------
# 3. Read BED files
# -----------------------------
mouse_specific <- import(mouse_specific_file)
human_specific <- import(human_specific_file)
conserved_hm   <- import(conserved_hm_file)

mouse_universe <- import(mouse_bg_file)
human_universe <- import(human_bg_file)

# -----------------------------
# 4. Clean helper
# -----------------------------
clean_gr <- function(gr) {
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)
  gr <- reduce(gr)
  gr
}

mouse_specific <- clean_gr(mouse_specific)
human_specific <- clean_gr(human_specific)
conserved_hm   <- clean_gr(conserved_hm)
mouse_universe <- clean_gr(mouse_universe)
human_universe <- clean_gr(human_universe)

# -----------------------------
# 5. Basic QC
# -----------------------------
mouse_specific_in_bg <- sum(countOverlaps(mouse_specific, mouse_universe) > 0)
human_specific_in_bg <- sum(countOverlaps(human_specific, human_universe) > 0)
conserved_hm_in_bg   <- sum(countOverlaps(conserved_hm, mouse_universe) > 0)

cat("=== Region counts after reduce() ===\n")
cat("mouse_specific:", length(mouse_specific), "\n")
cat("human_specific:", length(human_specific), "\n")
cat("conserved_human_in_mouse:", length(conserved_hm), "\n")
cat("mouse_universe:", length(mouse_universe), "\n")
cat("human_universe:", length(human_universe), "\n\n")

cat("=== Total widths ===\n")
cat("mouse_specific width:", sum(width(mouse_specific)), "\n")
cat("human_specific width:", sum(width(human_specific)), "\n")
cat("conserved_human_in_mouse width:", sum(width(conserved_hm)), "\n")
cat("mouse_universe width:", sum(width(mouse_universe)), "\n")
cat("human_universe width:", sum(width(human_universe)), "\n\n")

cat("=== Overlap with background ===\n")
cat("mouse_specific overlapping mouse_universe:", mouse_specific_in_bg, "of", length(mouse_specific), "\n")
cat("human_specific overlapping human_universe:", human_specific_in_bg, "of", length(human_specific), "\n")
cat("conserved_human_in_mouse overlapping mouse_universe:", conserved_hm_in_bg, "of", length(conserved_hm), "\n\n")

# -----------------------------
# 6. Run rGREAT
# -----------------------------
mouse_job <- great(
  gr = mouse_specific,
  gene_sets = ontology_to_run,
  tss_source = mouse_tss_source,
  background = mouse_universe
)

human_job <- great(
  gr = human_specific,
  gene_sets = ontology_to_run,
  tss_source = human_tss_source,
  background = human_universe
)

conserved_hm_job <- great(
  gr = conserved_hm,
  gene_sets = ontology_to_run,
  tss_source = mouse_tss_source,
  background = mouse_universe
)

# -----------------------------
# 7. Extract enrichment tables
# -----------------------------
mouse_tbl <- getEnrichmentTable(mouse_job)
human_tbl <- getEnrichmentTable(human_job)
conserved_hm_tbl <- getEnrichmentTable(conserved_hm_job)

write.csv(
  mouse_tbl,
  file.path(outdir, "mouse_specific_rGREAT_GO_BP.csv"),
  row.names = FALSE
)

write.csv(
  human_tbl,
  file.path(outdir, "human_specific_rGREAT_GO_BP.csv"),
  row.names = FALSE
)

write.csv(
  conserved_hm_tbl,
  file.path(outdir, "conserved_human_in_mouse_rGREAT_GO_BP.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Standardize result tables
# -----------------------------
standardize_tbl <- function(tb, label) {
  tb2 <- tb
  
  term_col <- intersect(c("name", "term_name", "description"), colnames(tb2))
  if (length(term_col) == 0) stop("Cannot find term name column.")
  term_col <- term_col[1]
  
  padj_col <- intersect(c("p_adjust", "adj_p_value", "fdr", "q_value"), colnames(tb2))
  if (length(padj_col) == 0) stop("Cannot find adjusted p-value column.")
  padj_col <- padj_col[1]
  
  enrich_col <- intersect(c("fold_enrichment", "fold_enrichment_binom", "enrichment"), colnames(tb2))
  if (length(enrich_col) == 0) {
    tb2$fold_enrichment_final <- NA_real_
  } else {
    enrich_col <- enrich_col[1]
    tb2$fold_enrichment_final <- tb2[[enrich_col]]
  }
  
  tb2$term_name_final <- tb2[[term_col]]
  tb2$padj_final <- tb2[[padj_col]]
  tb2$neglog10_padj <- -log10(tb2$padj_final)
  tb2$dataset <- label
  
  tb2
}

mouse_tbl2 <- standardize_tbl(mouse_tbl, "mouse_specific")
human_tbl2 <- standardize_tbl(human_tbl, "human_specific")
conserved_hm_tbl2 <- standardize_tbl(conserved_hm_tbl, "conserved_human_in_mouse")

# -----------------------------
# 9. Significant terms
# -----------------------------
mouse_sig <- mouse_tbl2 %>%
  filter(!is.na(padj_final), padj_final < 0.05) %>%
  arrange(padj_final)

human_sig <- human_tbl2 %>%
  filter(!is.na(padj_final), padj_final < 0.05) %>%
  arrange(padj_final)

conserved_hm_sig <- conserved_hm_tbl2 %>%
  filter(!is.na(padj_final), padj_final < 0.05) %>%
  arrange(padj_final)

write.csv(
  mouse_sig,
  file.path(outdir, "mouse_specific_rGREAT_GO_BP_sig.csv"),
  row.names = FALSE
)

write.csv(
  human_sig,
  file.path(outdir, "human_specific_rGREAT_GO_BP_sig.csv"),
  row.names = FALSE
)

write.csv(
  conserved_hm_sig,
  file.path(outdir, "conserved_human_in_mouse_rGREAT_GO_BP_sig.csv"),
  row.names = FALSE
)

cat("=== Significant term counts (FDR < 0.05) ===\n")
cat("mouse_specific:", nrow(mouse_sig), "\n")
cat("human_specific:", nrow(human_sig), "\n")
cat("conserved_human_in_mouse:", nrow(conserved_hm_sig), "\n\n")

# -----------------------------
# 10. Top 10 terms for plotting
# -----------------------------
mouse_plot <- mouse_tbl2 %>%
  filter(!is.na(padj_final), is.finite(neglog10_padj)) %>%
  arrange(padj_final) %>%
  slice_head(n = 10)

human_plot <- human_tbl2 %>%
  filter(!is.na(padj_final), is.finite(neglog10_padj)) %>%
  arrange(padj_final) %>%
  slice_head(n = 10)

conserved_hm_plot <- conserved_hm_tbl2 %>%
  filter(!is.na(padj_final), is.finite(neglog10_padj)) %>%
  arrange(padj_final) %>%
  slice_head(n = 10)

write.csv(
  mouse_plot,
  file.path(outdir, "mouse_specific_top10_for_plot.csv"),
  row.names = FALSE
)

write.csv(
  human_plot,
  file.path(outdir, "human_specific_top10_for_plot.csv"),
  row.names = FALSE
)

write.csv(
  conserved_hm_plot,
  file.path(outdir, "conserved_human_in_mouse_top10_for_plot.csv"),
  row.names = FALSE
)

# -----------------------------
# 11. Single-panel dotplot helper
# -----------------------------
make_single_dotplot <- function(tb, panel_title, filename, width = 7.5, height = 5.8) {
  if (nrow(tb) == 0) {
    message("No rows available for plotting: ", filename)
    return(invisible(NULL))
  }
  
  tb2 <- tb %>%
    arrange(padj_final) %>%
    mutate(
      term_name_final = str_wrap(term_name_final, width = 32),
      term_plot = factor(term_name_final, levels = rev(term_name_final))
    )
  
  max_x <- max(tb2$neglog10_padj, na.rm = TRUE)
  x_upper <- ceiling(max_x * 1.08 * 10) / 10
  x_upper <- max(x_upper, 1)
  
  p <- ggplot(tb2, aes(x = neglog10_padj, y = term_plot)) +
    geom_point(
      aes(size = fold_enrichment_final, fill = neglog10_padj),
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_x_continuous(limits = c(0, x_upper)) +
    scale_fill_gradient(
      low = "#DCE6F2",
      high = "#4C78A8",
      name = expression(-log[10]("adj. p"))
    ) +
    scale_size_continuous(
      range = c(3, 9),
      name = "Fold enrichment"
    ) +
    labs(
      title = panel_title,
      x = expression(-log[10]("adjusted p-value")),
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.y        = element_text(size = 10.5, color = "black"),
      axis.text.x        = element_text(size = 10.5, color = "black"),
      axis.title.x       = element_text(size = 11.5, color = "black"),
      axis.line.x        = element_line(color = "black", linewidth = 0.4),
      axis.line.y        = element_line(color = "black", linewidth = 0.4),
      axis.ticks         = element_line(color = "black", linewidth = 0.35),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.title       = element_text(size = 10.5),
      legend.text        = element_text(size = 9.5),
      legend.position    = "right",
      plot.margin        = margin(8, 16, 8, 8)
    )
  
  ggsave(
    file.path(outdir, filename),
    p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

# -----------------------------
# 12. Plot three separate figures
# -----------------------------
make_single_dotplot(
  tb = mouse_plot,
  panel_title = "mouse_specific",
  filename = "mouse_specific_top10_dotplot.png"
)

make_single_dotplot(
  tb = human_plot,
  panel_title = "human_specific",
  filename = "human_specific_top10_dotplot.png"
)

make_single_dotplot(
  tb = conserved_hm_plot,
  panel_title = "conserved_human_in_mouse",
  filename = "conserved_human_in_mouse_top10_dotplot.png"
)

# -----------------------------
# 13. Save combined significant table
# -----------------------------
combined_sig <- bind_rows(mouse_sig, human_sig, conserved_hm_sig)
write.csv(
  combined_sig,
  file.path(outdir, "combined_mouse_human_conserved_sig_terms.csv"),
  row.names = FALSE
)

cat("Done. Results saved in:", outdir, "\n")
