# =========================================================
# Task 5: Promoter/enhancer classification for OCR categories
# Compare species-specific and conserved peak sets using TSS-centered windows
# =========================================================

# -----------------------------
# 0. Load packages
# -----------------------------
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyr)

# -----------------------------
# 1. User settings and output paths
# -----------------------------
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
input_dir <- Sys.getenv("TASK5_INPUT_DIR", unset = file.path(project_root, "results", "mapping"))
outdir <- Sys.getenv("TASK5_OUTDIR", unset = file.path(project_root, "results", "task_5_enhancer_promoter"))
# Use environment variables so the same script can run from different working
# directories or be reused in a cluster job without editing file paths.

mouse_peak_file <- file.path(input_dir, "mouse_specific.bed")
human_peak_file <- file.path(input_dir, "human_specific.bed")
cons_human_in_mouse_file <- file.path(input_dir, "conserved_human_in_mouse.bed")
human_tss_file <- file.path(input_dir, "human_tss.bed")
mouse_tss_file <- file.path(input_dir, "mouse_tss.bed")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Check required input files
# -----------------------------
required_files <- c(
  mouse_peak_file,
  human_peak_file,
  cons_human_in_mouse_file,
  human_tss_file,
  mouse_tss_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required input files for Task 5:\n",
      paste(missing_files, collapse = "\n"),
      "\nSet TASK5_INPUT_DIR or PROJECT_ROOT as needed."
    )
  )
}

# -----------------------------
# 3. Helper functions for BED/TSS input
# -----------------------------
read_clean_bed <- function(file) {
  # Import peak intervals, drop non-standard chromosomes, and sort them so
  # overlap checks are based on a clean and consistent interval set.
  gr <- import(file)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)
  gr
}

read_teacher_tss <- function(file) {
  tb <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  # Convert the BED-like TSS table into GRanges. The start position is shifted
  # by +1 so the interval is interpreted correctly in 1-based GRanges space.
  gr <- GRanges(
    seqnames = tb[[1]],
    ranges = IRanges(start = tb[[2]] + 1, end = tb[[3]]),
    strand = tb[[5]]
  )
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)
  gr
}

# -----------------------------
# 4. Read peak sets and TSS annotations
# -----------------------------
mouse_specific <- read_clean_bed(mouse_peak_file)
human_specific <- read_clean_bed(human_peak_file)
conserved_human_in_mouse <- read_clean_bed(cons_human_in_mouse_file)

cat("mouse_specific:", length(mouse_specific), "\n")
cat("human_specific:", length(human_specific), "\n")
cat("conserved_human_in_mouse:", length(conserved_human_in_mouse), "\n\n")

cat("mouse_specific reduced:", length(reduce(mouse_specific)), "\n")
cat("human_specific reduced:", length(reduce(human_specific)), "\n")
cat("conserved_human_in_mouse reduced:", length(reduce(conserved_human_in_mouse)), "\n\n")

human_tss <- read_teacher_tss(human_tss_file)
mouse_tss <- read_teacher_tss(mouse_tss_file)

# -----------------------------
# 5. Build promoter windows
# -----------------------------
# promoter window = TSS +/- 2 kb
# Each TSS is expanded into a promoter-centered window so peaks can be assigned
# using a simple overlap test instead of nearest-distance logic.
human_promoters <- promoters(human_tss, upstream = 2000, downstream = 2000)
mouse_promoters <- promoters(mouse_tss, upstream = 2000, downstream = 2000)

human_promoters <- trim(human_promoters)
mouse_promoters <- trim(mouse_promoters)

human_promoters <- reduce(human_promoters)
mouse_promoters <- reduce(mouse_promoters)

# -----------------------------
# 6. Classify peaks as promoter or enhancer
# promoter = overlaps promoter window
# enhancer = does not overlap promoter window
# note: "Enhancer" here means non-promoter / distal peak by operational definition
# -----------------------------
classify_peak_set <- function(peaks, promoter_windows, dataset_name) {
  is_promoter <- overlapsAny(peaks, promoter_windows)

  # Keep one row per peak with a simple binary label used for summaries
  # and downstream plotting.
  tibble(
    dataset = dataset_name,
    seqnames = as.character(seqnames(peaks)),
    start = start(peaks),
    end = end(peaks),
    width = width(peaks),
    region_class = ifelse(is_promoter, "Promoter", "Enhancer")
  )
}

df_mouse_specific <- classify_peak_set(
  mouse_specific, mouse_promoters, "Mouse-specific"
)

df_human_specific <- classify_peak_set(
  human_specific, human_promoters, "Human-specific"
)

df_cons_hm <- classify_peak_set(
  conserved_human_in_mouse, mouse_promoters, "Conserved (human in mouse)"
)

df <- bind_rows(df_mouse_specific, df_human_specific, df_cons_hm)

# Fix the display order so summary tables and plots stay consistent across runs.
df$dataset <- factor(
  df$dataset,
  levels = c(
    "Mouse-specific",
    "Human-specific",
    "Conserved (human in mouse)"
  )
)

df$region_class <- factor(df$region_class, levels = c("Promoter", "Enhancer"))

# -----------------------------
# 7. Build summary tables
# -----------------------------
# Count how many peaks from each dataset fall into promoter vs enhancer bins.
count_df <- df %>%
  count(dataset, region_class, name = "count")

prop_df <- count_df %>%
  group_by(dataset) %>%
  mutate(
    total = sum(count),
    prop = count / total,
    prop_label = percent(prop, accuracy = 1)
  ) %>%
  # Put enhancer first so it is drawn at the bottom of the stacked bars and
  # promoter sits on top in a stable order.
  arrange(dataset, desc(region_class)) %>%   # Enhancer先（底部），Promoter后（顶部）
  group_by(dataset) %>%
  mutate(
    # Precompute label positions at the center of each stacked segment.
    label_y = cumsum(prop) - prop / 2
  ) %>%
  ungroup()

write.csv(df, file.path(outdir, "task5_peak_assignment.csv"), row.names = FALSE)
write.csv(count_df, file.path(outdir, "task5_count_summary.csv"), row.names = FALSE)
write.csv(prop_df, file.path(outdir, "task5_proportion_summary.csv"), row.names = FALSE)

# -----------------------------
# 8. Plot setup
# -----------------------------
plotdir <- file.path(outdir, "plots")
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

theme_clean <- function(base_size = 12) {
  # Shared theme keeps the three panels visually aligned while leaving
  # per-plot scales and legends free to differ when needed.
  theme_classic(base_size = base_size) +
    theme(
      axis.text    = element_text(color = "black"),
      axis.title   = element_text(color = "black"),
      axis.line    = element_line(linewidth = 0.45, color = "black"),
      axis.ticks   = element_line(linewidth = 0.4, color = "black"),
      legend.title = element_blank(),
      legend.text  = element_text(color = "black"),
      plot.title   = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      panel.border = element_blank()
    )
}

region_colors <- c("Promoter" = "#4C78A8", "Enhancer" = "#D98E73")

# -----------------------------
# 9. Plot 1: proportion stacked bar
# -----------------------------
# This panel emphasizes composition within each dataset rather than raw counts.
p1 <- ggplot(prop_df, aes(x = dataset, y = prop, fill = region_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = prop_label), size = 3.8, color = "black") +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0),
    limits = c(0, 1.01)
  ) +
  labs(title = "Peak composition", x = NULL, y = "Proportion of peaks") +
  theme_clean(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "top")

# -----------------------------
# 10. Plot 2: count grouped bar
# -----------------------------
# This panel preserves absolute counts so differences in dataset size remain visible.
p2 <- ggplot(count_df, aes(x = dataset, y = count, fill = region_class)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = comma(count)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(title = "Promoter vs enhancer counts", x = NULL, y = "Number of peaks") +
  theme_clean(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "top")

# -----------------------------
# 11. Plot 3: enhancer:promoter ratio
# -----------------------------
ratio_df <- count_df %>%
  pivot_wider(names_from = region_class, values_from = count) %>%
  mutate(
    # Replace missing classes with zero so ratios can still be computed when
    # a dataset contains only one category after classification.
    Promoter = ifelse(is.na(Promoter), 0, Promoter),
    Enhancer = ifelse(is.na(Enhancer), 0, Enhancer),
    # Avoid division by zero by returning NA when no promoter peaks are present.
    enhancer_promoter_ratio = ifelse(Promoter > 0, Enhancer / Promoter, NA_real_)
  )

p3 <- ggplot(ratio_df, aes(x = dataset, y = enhancer_promoter_ratio, fill = dataset)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(enhancer_promoter_ratio), "NA", round(enhancer_promoter_ratio, 1))),
    vjust = -0.4,
    size = 3.8,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "Mouse-specific" = "#7FB3D3",
    "Human-specific" = "#4C78A8",
    "Conserved (human in mouse)" = "#E8A87C"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Enhancer : Promoter ratio", x = NULL, y = "Enhancer / Promoter") +
  theme_clean(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "none")

# -----------------------------
# 12. Save final figure
# -----------------------------
# Save only the assembled patchwork figure; the component plots remain as
# in-memory objects used to build the final panel.
ggsave(
  file.path(plotdir, "all_three_panel.png"),
  p1 | p2 | p3,
  width = 18,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

cat("\nDone. Plots saved to:", plotdir, "\n")
