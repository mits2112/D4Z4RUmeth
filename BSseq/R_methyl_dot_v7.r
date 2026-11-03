#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
})

# ============================================================
# Usage:
#   Rscript R_methyl_dot_v7.r ${INDIR} ${OUTDIR} ${custom_fasta} ${condition_list} [facet_ncol]
#
# Arguments:
#   INDIR         Directory containing Bismark *.cov or *.cov.gz files
#   OUTDIR        Directory where output tables and figures are saved
#   custom_fasta  Custom reference FASTA used for Bismark mapping
#   condition_list Tab-separated file with columns: full_name, short_name
#                 full_name is used to find ${INDIR}/${full_name}.local_roi.cov.gz
#                 short_name is used as the plot label.
#   facet_ncol    Optional. Number of facet columns for
#                 region_methylation_dotplot_selected_samples_ordered_horizontal.png
#                 Default: 5
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript R_methyl_dot_v7.r <cov_dir> <out_dir> <custom_fasta> <condition_list.tsv> [facet_ncol]\n",
    "Example: Rscript R_methyl_dot_v7.r ${INDIR} ${OUTDIR} ${custom_fasta} ${OUTDIR}/condition_list.tsv 5"
  )
}

cov_dir <- args[1]
out_dir <- args[2]
fasta_file <- args[3]

condition_file <- args[4]
facet_ncol <- 5

if (length(args) >= 5) {
  if (!str_detect(args[5], "^[0-9]+$")) {
    stop("facet_ncol must be a positive integer.")
  }
  facet_ncol <- as.integer(args[5])
}

if (length(args) > 5) {
  stop(
    "Too many arguments. Usage: Rscript R_methyl_dot_v7.r ",
    "<cov_dir> <out_dir> <custom_fasta> <condition_list.tsv> [facet_ncol]"
  )
}

if (is.na(facet_ncol) || facet_ncol < 1) {
  stop("facet_ncol must be a positive integer.")
}

if (!file.exists(condition_file)) {
  stop("condition_list file does not exist: ", condition_file)
}

# The condition list explicitly defines which samples/files are used and how they are labeled.
# Required columns:
#   full_name   used to find cov file: <cov_dir>/<full_name>.local_roi.cov.gz
#   short_name  label used in figures
#
# Headerless two-column TSV is also accepted and interpreted as:
#   full_name <tab> short_name
read_condition_list <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!str_detect(lines, "^\\s*$")]
  lines <- lines[!str_detect(lines, "^\\s*#")]

  if (length(lines) == 0) {
    stop("condition_list is empty: ", path)
  }

  first_fields <- str_split(lines[1], "\\t")[[1]]

  if (all(c("full_name", "short_name") %in% first_fields)) {
    x <- read_tsv(
      path,
      col_types = cols(.default = col_character()),
      progress = FALSE,
      comment = "#"
    ) |>
      transmute(
        full_name = str_trim(full_name),
        short_name = str_trim(short_name)
      )
  } else if (all(c("sample", "condition") %in% first_fields)) {
    # Backward-compatible import, if an old condition_list-style file is supplied.
    x <- read_tsv(
      path,
      col_types = cols(.default = col_character()),
      progress = FALSE,
      comment = "#"
    ) |>
      transmute(
        full_name = str_trim(sample),
        short_name = str_trim(condition)
      )
  } else {
    x <- read_tsv(
      path,
      col_names = c("full_name", "short_name"),
      col_types = cols(.default = col_character()),
      progress = FALSE,
      comment = "#"
    ) |>
      transmute(
        full_name = str_trim(full_name),
        short_name = str_trim(short_name)
      )
  }

  x |>
    filter(!is.na(full_name), !is.na(short_name), full_name != "", short_name != "")
}

condition_list <- read_condition_list(condition_file)

if (nrow(condition_list) == 0) {
  stop("condition_list has no valid rows after removing empty rows.")
}

duplicated_full_names <- condition_list |>
  count(full_name, name = "n") |>
  filter(n > 1) |>
  pull(full_name)

if (length(duplicated_full_names) > 0) {
  stop(
    "Each full_name must appear only once in condition_list. Duplicated full_name(s): ",
    paste(duplicated_full_names, collapse = ", ")
  )
}

# Plot order follows the row order in condition_list.tsv.
samples <- condition_list$full_name
sample_order <- condition_list$short_name
sample_order_unique <- unique(sample_order)

if (length(sample_order_unique) < length(sample_order)) {
  warning(
    "Some short_name values are duplicated in condition_list.tsv. ",
    "The selected-sample grouped plot will merge rows with the same short_name."
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(cov_dir)) {
  stop("Input cov_dir does not exist: ", cov_dir)
}
if (!file.exists(fasta_file)) {
  stop("custom_fasta does not exist: ", fasta_file)
}

# -------------------------
# parameters
# -------------------------
min_cov <- 1
keep_fraction <- 0.51  # 0.5
max_gap <- 100
min_cpg_per_region <- 5  # 3
min_region_cov <- 10  # 5

message("Input cov_dir: ", cov_dir)
message("Output out_dir: ", out_dir)
message("Custom FASTA: ", fasta_file)
message("Condition list: ", condition_file)
message("Selected full names from condition_list: ", paste(samples, collapse = ", "))
message("Plot labels from condition_list: ", paste(sample_order, collapse = ", "))
message("Facet ncol for horizontal selected plot: ", facet_ncol)

# -------------------------
# 1. Extract all CpGs from custom reference FASTA
# -------------------------
ref <- readDNAStringSet(fasta_file)
ref_lengths <- width(ref)
custom_genome_length <- max(ref_lengths)

all_cpg <- imap_dfr(as.character(ref), function(seq, chr) {
  pos <- gregexpr("CG", seq, fixed = TRUE)[[1]]

  if (pos[1] == -1) {
    return(tibble())
  }

  tibble(
    chr = chr,
    cpg_pos = as.integer(pos),
    locus = paste(chr, cpg_pos, sep = ":")
  )
})

message("Number of CpGs in reference: ", nrow(all_cpg))

if (nrow(all_cpg) == 0) {
  stop("No CpG sites found in custom_fasta: ", fasta_file)
}

# -------------------------
# 1b. Plot all CpG sites on the custom reference coordinates
# -------------------------
ref_coordinate_df <- tibble(
  chr = names(ref),
  seqlen = as.integer(width(ref)),
  y = seq_along(names(ref))
)

all_cpg_plot_df <- all_cpg |>
  left_join(ref_coordinate_df |> select(chr, y), by = "chr")

p_all_cpg <- ggplot() +
  geom_segment(
    data = ref_coordinate_df,
    aes(x = 1, xend = seqlen, y = y, yend = y),
    linewidth = 0.35
  ) +
  geom_segment(
    data = all_cpg_plot_df,
    aes(x = cpg_pos, xend = cpg_pos, y = y - 0.18, yend = y + 0.18),
    linewidth = 0.25,
    alpha = 0.75
  ) +
  scale_x_continuous(
    limits = c(1, max(ref_coordinate_df$seqlen)),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = ref_coordinate_df$y,
    labels = ref_coordinate_df$chr
  ) +
  labs(
    x = "Position on custom genome (bp)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "all_CpG_sites_on_custom_genome.png"),
  plot = p_all_cpg,
  width = 10,
  height = max(2.5, nrow(ref_coordinate_df) * 0.45 + 1.5),
  dpi = 300
)

# -------------------------
# 2. Read selected Bismark cov files defined in condition_list.tsv
# -------------------------
condition_file_table <- condition_list |>
  mutate(
    sample = full_name,
    sample_group = short_name,
    expected_cov_file = file.path(cov_dir, paste0(full_name, ".local_roi.cov.gz")),
    cov_file = expected_cov_file,
    in_condition_list = TRUE,
    in_cov_files = file.exists(cov_file),
    match_status = if_else(in_cov_files, "assigned", "cov_file_not_found"),
    match_rule = "condition_list_exact_file_name"
  )

sample_condition_assignment <- condition_file_table |>
  select(
    sample,
    sample_group,
    full_name,
    short_name,
    cov_file,
    expected_cov_file,
    in_condition_list,
    in_cov_files,
    match_status,
    match_rule
  ) |>
  mutate(
    sample_group = factor(sample_group, levels = sample_order_unique),
    sample_group = as.character(sample_group)
  ) |>
  arrange(match_status != "assigned", match(full_name, samples))

write_tsv(
  sample_condition_assignment,
  file.path(out_dir, "sample_condition_assignment.tsv")
)

sample_condition_assignment_summary <- sample_condition_assignment |>
  count(sample_group, match_status, match_rule, in_condition_list, in_cov_files, name = "n_samples") |>
  arrange(is.na(sample_group), sample_group, match_status, match_rule)

write_tsv(
  sample_condition_assignment_summary,
  file.path(out_dir, "sample_condition_assignment_summary.tsv")
)

missing_cov_files <- sample_condition_assignment |>
  filter(!in_cov_files) |>
  pull(expected_cov_file)

if (length(missing_cov_files) > 0) {
  stop(
    "Some cov files listed by condition_list.tsv were not found. ",
    "Expected file pattern is <INDIR>/<full_name>.local_roi.cov.gz. ",
    "See sample_condition_assignment.tsv. Missing file(s): ",
    paste(missing_cov_files, collapse = ", ")
  )
}

cov_files <- condition_file_table$cov_file
names(cov_files) <- condition_file_table$full_name

message("Number of selected cov files: ", length(cov_files))

read_bismark_cov <- function(file, sample_name) {
  read_tsv(
    file,
    col_names = c("chr", "start", "end", "meth_pct", "mC", "unmC"),
    col_types = cols(
      chr = col_character(),
      start = col_integer(),
      end = col_integer(),
      meth_pct = col_double(),
      mC = col_double(),
      unmC = col_double()
    ),
    progress = FALSE
  ) |>
    transmute(
      sample = sample_name,
      chr,
      start,
      end,
      cpg_pos = start,
      locus = paste(chr, start, sep = ":"),
      locus_with_end = paste(chr, start, end, sep = ":"),
      meth_pct,
      mC,
      unmC,
      coverage = mC + unmC,
      cpg_meth_pct = if_else(coverage > 0, 100 * mC / coverage, NA_real_)
    )
}

cov_df <- map2_dfr(cov_files, names(cov_files), read_bismark_cov)

if (nrow(cov_df) == 0) {
  stop("No rows were read from selected cov files.")
}

# Keep the sample order exactly as specified in condition_list.tsv.
samples <- condition_file_table$full_name
message("Number of selected samples: ", length(samples))

# -------------------------
# 3. Make all-CpG x all-sample coverage table
#    CpGs absent from Bismark cov output are treated as coverage 0.
# -------------------------
coverage_table <- all_cpg |>
  crossing(sample = samples) |>
  left_join(
    cov_df |>
      select(sample, chr, cpg_pos, locus, meth_pct, mC, unmC, coverage),
    by = c("chr", "cpg_pos", "locus", "sample")
  ) |>
  mutate(
    coverage = replace_na(coverage, 0),
    mC = replace_na(mC, 0),
    unmC = replace_na(unmC, 0),
    meth_pct = if_else(coverage > 0, meth_pct, NA_real_),
    pass_cov = coverage >= min_cov
  )

write_tsv(
  coverage_table,
  file.path(out_dir, "all_CpG_coverage_by_sample.tsv")
)

# -------------------------
# 4. Sample-level coverage summary
# -------------------------
sample_summary <- coverage_table |>
  group_by(sample) |>
  summarise(
    n_total_cpg = n(),
    n_covered_cpg = sum(coverage > 0),
    n_pass_cpg = sum(coverage >= min_cov),
    mean_cov = mean(coverage),
    median_cov = median(coverage),
    .groups = "drop"
  )

write_tsv(
  sample_summary,
  file.path(out_dir, "sample_coverage_summary.tsv")
)

# -------------------------
# 5. Coverage threshold summary
# -------------------------
thresholds <- c(1, 2, 3, 4, 5, 10, 20, 30, 50)

threshold_summary <- map_dfr(thresholds, function(th) {
  coverage_table |>
    group_by(sample) |>
    summarise(
      min_cov = th,
      n_pass_cpg = sum(coverage >= th),
      fraction_pass_cpg = n_pass_cpg / n(),
      .groups = "drop"
    )
})

write_tsv(
  threshold_summary,
  file.path(out_dir, "coverage_threshold_summary_by_sample.tsv")
)

p <- ggplot(threshold_summary, aes(x = min_cov, y = n_pass_cpg, color = sample)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Minimum coverage",
    y = "Number of usable CpGs",
    color = "Sample"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "coverage_threshold_summary_by_sample.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------
# 6. CpG-level coverage summary
# -------------------------
cpg_summary <- coverage_table |>
  group_by(chr, cpg_pos, locus) |>
  summarise(
    n_samples = n(),
    n_samples_pass = sum(coverage >= min_cov),
    fraction_samples_pass = n_samples_pass / n_samples,
    mean_cov = mean(coverage),
    median_cov = median(coverage),
    .groups = "drop"
  )

write_tsv(
  cpg_summary,
  file.path(out_dir, "CpG_coverage_summary.tsv")
)

usable_cpg_80 <- cpg_summary |>
  filter(fraction_samples_pass >= 0.8)

write_tsv(
  usable_cpg_80,
  file.path(out_dir, "usable_CpGs_minCov1_80percentSamples.tsv")
)

# -------------------------
# 7. CpG position coverage plots
# -------------------------
p1 <- ggplot(coverage_table, aes(x = cpg_pos, y = coverage, color = sample)) +
  geom_line(alpha = 0.8) +
  geom_hline(yintercept = min_cov, linetype = "dashed") +
  labs(
    x = "CpG position on custom genome",
    y = "Coverage",
    color = "Sample"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "CpG_position_coverage_lineplot.png"),
  plot = p1,
  width = 8,
  height = 5,
  dpi = 300
)

p2 <- ggplot(coverage_table, aes(x = cpg_pos, y = sample, fill = coverage)) +
  geom_tile() +
  geom_point(
    data = coverage_table |> filter(coverage >= min_cov),
    aes(x = cpg_pos, y = sample),
    inherit.aes = FALSE,
    size = 0.8
  ) +
  labs(
    x = "CpG position on custom genome",
    y = "Sample",
    fill = "Coverage"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "CpG_position_sample_coverage_heatmap.png"),
  plot = p2,
  width = 9,
  height = max(5, length(samples) * 0.3),
  dpi = 300
)

# -------------------------
# 8. Estimate usable regions from CpGs covered in >= keep_fraction samples
# -------------------------
usable_cpg <- cpg_summary |>
  filter(fraction_samples_pass >= keep_fraction) |>
  arrange(chr, cpg_pos)

if (nrow(usable_cpg) == 0) {
  stop("No usable CpGs remained. Consider lowering min_cov or keep_fraction.")
}

usable_regions <- usable_cpg |>
  arrange(chr, cpg_pos) |>
  group_by(chr) |>
  mutate(
    gap = cpg_pos - lag(cpg_pos),
    new_region = if_else(is.na(gap) | gap > max_gap, 1, 0),
    region_index = cumsum(new_region)
  ) |>
  group_by(chr, region_index) |>
  summarise(
    region_start = min(cpg_pos),
    region_end = max(cpg_pos) + 1,
    n_cpg = n(),
    mean_fraction_samples_pass = mean(fraction_samples_pass),
    mean_cov = mean(mean_cov),
    .groups = "drop"
  ) |>
  mutate(
    region_id = paste0(chr, ":", region_start, "-", region_end)
  ) |>
  select(
    region_id, chr, region_start, region_end, n_cpg,
    mean_fraction_samples_pass, mean_cov
  )

write_tsv(
  usable_regions,
  file.path(out_dir, "estimated_usable_regions.tsv")
)

usable_regions_bed <- usable_regions |>
  transmute(
    chr = chr,
    start0 = region_start - 1,
    end = region_end,
    region_id = region_id
  )

write_tsv(
  usable_regions_bed,
  file.path(out_dir, "estimated_usable_regions.bed"),
  col_names = FALSE
)

p3 <- ggplot() +
  geom_segment(
    data = tibble(x = 1, xend = custom_genome_length, y = 0, yend = 0),
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 1
  ) +
  geom_rect(
    data = usable_regions,
    aes(
      xmin = region_start,
      xmax = region_end,
      ymin = -0.1,
      ymax = 0.1
    ),
    alpha = 0.6
  ) +
  geom_point(
    data = cpg_summary,
    aes(x = cpg_pos, y = 0.25, size = n_samples_pass),
    alpha = 0.6
  ) +
  scale_x_continuous(limits = c(1, custom_genome_length)) +
  labs(
    x = "Position on custom genome",
    y = NULL,
    size = paste0("Samples with coverage >= ", min_cov)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "estimated_usable_regions.png"),
  plot = p3,
  width = 10,
  height = 4,
  dpi = 300
)

# -------------------------
# 9. CpG-level methylation dot plot
# -------------------------
n_samples <- n_distinct(cov_df$sample)
min_samples <- ceiling(n_samples * keep_fraction)

keep_loci <- cov_df |>
  mutate(pass = coverage >= min_cov) |>
  group_by(locus_with_end) |>
  summarise(
    n_pass_samples = sum(pass),
    .groups = "drop"
  ) |>
  filter(n_pass_samples >= min_samples)

plot_cpg_df <- cov_df |>
  semi_join(keep_loci, by = "locus_with_end") |>
  filter(coverage >= min_cov)

p <- ggplot(plot_cpg_df, aes(x = sample, y = meth_pct, size = coverage)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.35) +
  scale_size_area(max_size = 6) +
  labs(
    x = "Sample",
    y = "Methylation (%)",
    size = "Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_dir, "CpG_methylation_dotplot_by_sample.png"),
  plot = p,
  width = max(7, length(unique(plot_cpg_df$sample)) * 0.5),
  height = 5,
  dpi = 300
)

# -------------------------
# 10. Region-level methylation summaries with CpG-level variation
# -------------------------
regions <- usable_regions_bed |>
  mutate(start = start0 + 1)

gr_regions <- GRanges(
  seqnames = regions$chr,
  ranges = IRanges(start = regions$start, end = regions$end),
  region_id = regions$region_id
)

read_region_bismark_cov <- function(file, sample_name) {
  cov <- read_tsv(
    file,
    col_names = c("chr", "start", "end", "meth_pct", "mC", "unmC"),
    col_types = cols(
      chr = col_character(),
      start = col_integer(),
      end = col_integer(),
      meth_pct = col_double(),
      mC = col_double(),
      unmC = col_double()
    ),
    progress = FALSE
  ) |>
    mutate(
      sample = sample_name,
      coverage = mC + unmC,
      cpg_meth_pct = if_else(coverage > 0, 100 * mC / coverage, NA_real_)
    )

  gr_cpg <- GRanges(
    seqnames = cov$chr,
    ranges = IRanges(start = cov$start, end = cov$end),
    mC = cov$mC,
    unmC = cov$unmC,
    coverage = cov$coverage,
    cpg_meth_pct = cov$cpg_meth_pct,
    sample = cov$sample
  )

  hits <- findOverlaps(gr_regions, gr_cpg)

  tibble(
    region_id = mcols(gr_regions)$region_id[queryHits(hits)],
    sample = mcols(gr_cpg)$sample[subjectHits(hits)],
    mC = mcols(gr_cpg)$mC[subjectHits(hits)],
    unmC = mcols(gr_cpg)$unmC[subjectHits(hits)],
    coverage = mcols(gr_cpg)$coverage[subjectHits(hits)],
    cpg_meth_pct = mcols(gr_cpg)$cpg_meth_pct[subjectHits(hits)]
  ) |>
    group_by(region_id, sample) |>
    summarise(
      n_cpg = n(),
      total_mC = sum(mC),
      total_unmC = sum(unmC),
      total_cov = sum(coverage),

      # Coverage-weighted methylation across the region.
      meth_pct = 100 * total_mC / total_cov,

      # CpG-to-CpG variability within the region and sample.
      cpg_meth_min = min(cpg_meth_pct, na.rm = TRUE),
      cpg_meth_max = max(cpg_meth_pct, na.rm = TRUE),
      cpg_meth_mean = mean(cpg_meth_pct, na.rm = TRUE),
      cpg_meth_sd = sd(cpg_meth_pct, na.rm = TRUE),
      .groups = "drop"
    )
}

region_df <- map2_dfr(cov_files, names(cov_files), read_region_bismark_cov)

if (nrow(region_df) == 0) {
  stop("No CpGs overlapped estimated usable regions.")
}

n_samples_region <- n_distinct(region_df$sample)
min_samples_region <- ceiling(n_samples_region * keep_fraction)

keep_regions <- region_df |>
  mutate(
    pass = total_cov >= min_region_cov & n_cpg >= min_cpg_per_region
  ) |>
  group_by(region_id) |>
  summarise(
    n_pass_samples = sum(pass),
    .groups = "drop"
  ) |>
  filter(n_pass_samples >= min_samples_region)

plot_df <- region_df |>
  semi_join(keep_regions, by = "region_id") |>
  filter(
    total_cov >= min_region_cov,
    n_cpg >= min_cpg_per_region
  )

write_tsv(
  region_df,
  file.path(out_dir, "region_methylation_summary_all_regions.tsv")
)
write_tsv(
  plot_df,
  file.path(out_dir, "region_methylation_summary_plot_regions.tsv")
)

# Region methylation dot plot by sample.
p <- ggplot(plot_df, aes(x = sample, y = meth_pct, size = total_cov)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.6) +
  scale_size_area(max_size = 8) +
  labs(
    x = "Sample",
    y = "Methylation (%)",
    size = "Region coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_dir, "region_methylation_dotplot_by_sample.png"),
  plot = p,
  width = max(7, length(unique(plot_df$sample)) * 0.5),
  height = 5,
  dpi = 300
)

# NOTE:
# The old output "region_methylation_dotplot_by_sample_with_CpG_variation.png"
# has intentionally been removed.

# -------------------------
# 11. Selected-sample region plot with CpG-level SD bars
# -------------------------
plot_df_ordered <- plot_df |>
  left_join(
    sample_condition_assignment |>
      select(sample, sample_group),
    by = "sample"
  ) |>
  filter(!is.na(sample_group))

if (nrow(plot_df_ordered) > 0) {
  sample_order_present <- sample_order_unique[sample_order_unique %in% unique(plot_df_ordered$sample_group)]

  region_order <- regions |>
    arrange(chr, start0, end) |>
    pull(region_id) |>
    unique()

  plot_df_ordered <- plot_df_ordered |>
    mutate(
      sample_group = factor(sample_group, levels = sample_order_present),
      region_id = factor(region_id, levels = region_order),
      total_cov_for_size = pmin(total_cov, 300)
    )

  # Visualize the regions actually used in
  # region_methylation_dotplot_selected_samples_ordered_horizontal.png.
  selected_region_ids <- plot_df_ordered |>
    pull(region_id) |>
    as.character() |>
    unique()

  selected_regions_for_horizontal_plot <- regions |>
    filter(region_id %in% selected_region_ids) |>
    arrange(chr, start0, end)

  if (nrow(selected_regions_for_horizontal_plot) > 0) {
    ref_plot_df <- tibble(
      chr = names(ref),
      seqlen = as.integer(width(ref)),
      y = seq_along(names(ref))
    ) |>
      filter(chr %in% unique(selected_regions_for_horizontal_plot$chr))

    selected_regions_for_horizontal_plot <- selected_regions_for_horizontal_plot |>
      left_join(ref_plot_df |> select(chr, y), by = "chr") |>
      mutate(
        start = start0 + 1
      )

    p_region_used <- ggplot() +
      geom_segment(
        data = ref_plot_df,
        aes(x = 1, xend = seqlen, y = y, yend = y),
        linewidth = 0.35
      ) +
      geom_rect(
        data = selected_regions_for_horizontal_plot,
        aes(
          xmin = start,
          xmax = end,
          ymin = y - 0.12,
          ymax = y + 0.12
        ),
        alpha = 0.7
      ) +
      scale_x_continuous(
        limits = c(1, max(ref_plot_df$seqlen)),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(
        breaks = ref_plot_df$y,
        labels = ref_plot_df$chr
      ) +
      labs(
        x = "Position on custom genome (bp)",
        y = NULL
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.ticks.y = element_blank()
      )

    ggsave(
      filename = file.path(out_dir, "regions_used_in_region_methylation_dotplot_selected_samples_ordered_horizontal.png"),
      plot = p_region_used,
      width = 10,
      height = max(2, nrow(ref_plot_df) * 1.2),
      dpi = 300
    )
  } else {
    warning("No regions were used in the selected horizontal plot. Region overview plot was not generated.")
  }

  p <- ggplot(plot_df_ordered, aes(x = sample_group, y = meth_pct)) +
    geom_errorbar(
      aes(
        ymin = pmax(cpg_meth_mean - cpg_meth_sd, 0),
        ymax = pmin(cpg_meth_mean + cpg_meth_sd, 100),
        color = sample_group
      ),
      width = 0.08,
      linewidth = 0.25,
      alpha = 0.6
    ) +
    geom_point(
      aes(
        size = total_cov_for_size,
        color = sample_group
      ),
      alpha = 0.8
    ) +
    scale_size_area(
      max_size = 4,
      breaks = c(50, 100, 200, 300),
      labels = c("50", "100", "200", "300+")
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    facet_wrap(~ region_id) +
    labs(
      x = "Sample",
      y = "Methylation (%)",
      size = "Coverage",
      color = "Sample"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(out_dir, "region_methylation_dotplot_selected_samples_ordered.png"),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )

  p <- ggplot(plot_df_ordered, aes(x = meth_pct, y = sample_group)) +
    geom_errorbar(
      aes(
        xmin = pmax(cpg_meth_mean - cpg_meth_sd, 0),
        xmax = pmin(cpg_meth_mean + cpg_meth_sd, 100),
        color = sample_group
      ),
      width = 0.08,
      linewidth = 0.25,
      alpha = 0.6,
      orientation = "y"
    ) +
    geom_point(
      aes(
        size = total_cov_for_size,
        color = sample_group
      ),
      alpha = 0.8
    ) +
    scale_size_area(
      max_size = 4,
      breaks = c(50, 100, 200, 300),
      labels = c("50", "100", "200", "300+")
    ) +
    scale_x_continuous(
      limits = c(-10, 105),
      breaks = seq(0, 100, 25)
    ) +
    scale_y_discrete(limits = rev(sample_order_present)) +
    facet_wrap(
      ~ region_id,
      ncol = facet_ncol,
      strip.position = "top",
      labeller = label_wrap_gen(width = 12)
    ) +
    labs(
      x = "Methylation (%)",
      y = "Sample",
      size = "Coverage",
      color = "Sample"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 7, lineheight = 0.85, margin = margin(t = 2, b = 2)),
      plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
    )

  ggsave(
    filename = file.path(out_dir, "region_methylation_dotplot_selected_samples_ordered_horizontal.png"),
    plot = p,
    width = 10,
    height = 4,
    dpi = 300
  )
} else {
  warning("No selected sample groups were found. Selected-sample plots were not generated.")
}

message("Done. Outputs saved to: ", out_dir)
