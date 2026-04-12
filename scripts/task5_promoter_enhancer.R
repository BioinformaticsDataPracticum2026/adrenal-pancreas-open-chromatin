.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/4.4", .libPaths()))

library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

workdir <- "/jet/home/wzhang37/step5_enhancer_promoter"

mouse_peak_file <- file.path(workdir, "mouse_specific.bed")
human_peak_file <- file.path(workdir, "human_specific.bed")
cons_human_in_mouse_file <- file.path(workdir, "conserved_human_in_mouse.bed")
cons_mouse_in_human_file <- file.path(workdir, "conserved_mouse_in_human.bed")

human_tss_file <- "/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed"
mouse_tss_file <- "/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed"

outdir <- file.path(workdir, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_clean_bed <- function(file) {
  gr <- import(file)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)
  gr
}

read_teacher_tss <- function(file) {
  tb <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = tb[[1]],
    ranges = IRanges(start = tb[[2]] + 1, end = tb[[3]]),
    strand = tb[[5]]
  )
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)
  gr
}

mouse_specific <- read_clean_bed(mouse_peak_file)
human_specific <- read_clean_bed(human_peak_file)
conserved_human_in_mouse <- read_clean_bed(cons_human_in_mouse_file)
conserved_mouse_in_human <- read_clean_bed(cons_mouse_in_human_file)

human_tss <- read_teacher_tss(human_tss_file)
mouse_tss <- read_teacher_tss(mouse_tss_file)

# promoter window = TSS ± 2 kb
human_promoters <- promoters(human_tss, upstream = 2000, downstream = 2000)
mouse_promoters <- promoters(mouse_tss, upstream = 2000, downstream = 2000)

human_promoters <- trim(human_promoters)
mouse_promoters <- trim(mouse_promoters)

human_promoters <- reduce(human_promoters)
mouse_promoters <- reduce(mouse_promoters)

# -----------------------------
# classify peaks
# promoter = overlaps promoter window
# enhancer = does not overlap promoter window
# -----------------------------
classify_peak_set <- function(peaks, promoter_windows, dataset_name) {
  is_promoter <- overlapsAny(peaks, promoter_windows)

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

df_cons_mh <- classify_peak_set(
  conserved_mouse_in_human, human_promoters, "Conserved (mouse in human)"
)

df <- bind_rows(df_mouse_specific, df_human_specific, df_cons_hm, df_cons_mh)

df$dataset <- factor(
  df$dataset,
  levels = c(
    "Mouse-specific",
    "Human-specific",
    "Conserved (human in mouse)",
    "Conserved (mouse in human)"
  )
)

df$region_class <- factor(df$region_class, levels = c("Promoter", "Enhancer"))

# -----------------------------
# summaries
# -----------------------------
count_df <- df %>%
  count(dataset, region_class, name = "count")

prop_df <- count_df %>%
  group_by(dataset) %>%
  mutate(
    total = sum(count),
    prop = count / total,
    prop_label = percent(prop, accuracy = 1)
  ) %>%
  arrange(dataset, region_class) %>%
  group_by(dataset) %>%
  mutate(
    label_y = cumsum(prop) - prop / 2
  ) %>%
  ungroup()

write.csv(df, file.path(outdir, "step5_peak_assignment.csv"), row.names = FALSE)
write.csv(count_df, file.path(outdir, "step5_count_summary.csv"), row.names = FALSE)
write.csv(prop_df, file.path(outdir, "step5_proportion_summary.csv"), row.names = FALSE)

# -----------------------------
# plot theme
# -----------------------------
theme_clean <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5)
    )
}

region_colors <- c(
  "Promoter" = "#4C78A8",
  "Enhancer" = "#D98E73"
)

# -----------------------------
# plot 1: proportions
# -----------------------------
p1 <- ggplot(prop_df, aes(x = dataset, y = prop, fill = region_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = prop_label), size = 4) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Peak composition",
    x = NULL,
    y = "Proportion of peaks"
  ) +
  theme_clean(13) +
  theme(
    axis.text.x = element_text(angle = 16, hjust = 1),
    legend.position = "top"
  )

# -----------------------------
# plot 2: counts
# -----------------------------
p2 <- ggplot(count_df, aes(x = dataset, y = count, fill = region_class)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 4
  ) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Promoter vs enhancer counts",
    x = NULL,
    y = "Number of peaks"
  ) +
  theme_clean(13) +
  theme(
    axis.text.x = element_text(angle = 16, hjust = 1),
    legend.position = "top"
  )

final_plot <- p1 | p2

ggsave(
  file.path(outdir, "step5_promoter_enhancer_two_panel.png"),
  final_plot,
  width = 13,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(outdir, "step5_promoter_enhancer_two_panel.pdf"),
  final_plot,
  width = 13,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

print(final_plot)
