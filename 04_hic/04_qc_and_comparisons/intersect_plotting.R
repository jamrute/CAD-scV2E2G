library(dplyr)
library(data.table)
library(HiCcompare)
library(multiHiCcompare)
library(SpectralTAD)
library(TADCompare)
library(BiocParallel)
library(ggplot2)
library(Hmisc)
library(stringr)
library(pheatmap)


setwd("<path>/FULL_SAMPLES/INTERSECTS/")

#Get all sample IDs in a vector
samples <- c(SAMPLE_05, SAMPLE_09, SAMPLE_10, SAMPLE_11, SAMPLE_12, SAMPLE_13, SAMPLE_14, SAMPLE_15, SAMPLE_16, SAMPLE_17, SAMPLE_18, SAMPLE_19, SAMPLE_20, SAMPLE_21, SAMPLE_22)
resolutions <- c(10, 5, 2)

overlaps <- data.frame()

for (sample in samples) {
  for (resolution in resolutions){
    tmp.df <- fread(paste0(sample, "_", resolution,"kb_loopcaller_intersects.txt"))
    colnames(tmp.df) <- c("FitHiC2", "Shared", "HiCCUPS")
    tmp.df$sample <- sample
    tmp.df$resolution <- resolution*1000
    tmp.df <- melt(tmp.df, id.vars=c("sample", "resolution"))
    overlaps <- rbind(overlaps, tmp.df)
  }
}

overlaps$sample <- as.factor(overlaps$sample)

#10kb plot
ggplot(data=dplyr::filter(overlaps, resolution==10000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Significant Loop Count") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 10kb") + guides(fill=guide_legend(title="Loop Origin"))
ggsave("../QCFIGS/intra_method_loops_10kb_counts.pdf")
ggplot(data=dplyr::filter(overlaps, resolution==10000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Percentage of Significant Loops") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 10kb") + guides(fill=guide_legend(title="Loop Origin")) + 
  scale_y_continuous(labels = scales::percent_format())
ggsave("../QCFIGS/intra_method_loops_10kb_percentages.pdf")

#5kb plot
ggplot(data=dplyr::filter(overlaps, resolution==5000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Significant Loop Count") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 5kb") + guides(fill=guide_legend(title="Loop Origin"))
ggsave("../QCFIGS/intra_method_loops_5kb_counts.pdf")
ggplot(data=dplyr::filter(overlaps, resolution==5000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Percentage of Significant Loops") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 5kb") + guides(fill=guide_legend(title="Loop Origin")) + 
  scale_y_continuous(labels = scales::percent_format())
ggsave("../QCFIGS/intra_method_loops_5kb_percentages.pdf")

#2kb plot
ggplot(data=dplyr::filter(overlaps, resolution==2000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Significant Loop Count") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 2kb") + guides(fill=guide_legend(title="Loop Origin"))
ggsave("../QCFIGS/intra_method_loops_2kb_counts.pdf")
ggplot(data=dplyr::filter(overlaps, resolution==2000), aes(x=sample, group=sample, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab("Sample ID") + ylab("Percentage of Significant Loops") + ggtitle("Intra-Sample Comparison of Significant Loops, FitHiC2 vs. HiCCUPS, 2kb") + guides(fill=guide_legend(title="Loop Origin")) + 
  scale_y_continuous(labels = scales::percent_format())
ggsave("../QCFIGS/intra_method_loops_2kb_percentages.pdf")


### Go ahead and also plot out contact distance distributions per-sample for FitHiC2 and HiCCUPS.
library(data.table)
library(ggplot2)
library(dplyr)

# Samples/resolutions from your earlier block
samples <- c(SAMPLE_05, SAMPLE_09, SAMPLE_10, SAMPLE_11, SAMPLE_12, SAMPLE_13, SAMPLE_14, SAMPLE_15,
             SAMPLE_16, SAMPLE_17, SAMPLE_18, SAMPLE_19, SAMPLE_20, SAMPLE_21, SAMPLE_22)
resolutions <- c(10000, 5000, 2000)

# Toggle display cap without removing data
cap_y_at_2e6 <- TRUE
y_cap <- 2e6

# Helper: read BEDPE-like loop file and compute contact distance
get_contact_distances <- function(file, sample_id, resolution, method, header = TRUE) {
  dt <- fread(file, header = header, sep = "\t", quote = "", data.table = FALSE)
  
  if (ncol(dt) < 6) {
    stop(paste("File has fewer than 6 columns:", file))
  }
  
  dt <- dt[, 1:6]
  colnames(dt) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  
  dt$chr1 <- gsub("^chr", "", dt$chr1)
  dt$chr2 <- gsub("^chr", "", dt$chr2)
  
  dt$start1 <- as.numeric(dt$start1)
  dt$end1   <- as.numeric(dt$end1)
  dt$start2 <- as.numeric(dt$start2)
  dt$end2   <- as.numeric(dt$end2)
  
  dt <- dt %>%
    filter(!is.na(start1), !is.na(end1), !is.na(start2), !is.na(end2)) %>%
    filter(chr1 == chr2)
  
  dt$mid1 <- (dt$start1 + dt$end1) / 2
  dt$mid2 <- (dt$start2 + dt$end2) / 2
  dt$contact_distance_bp <- abs(dt$mid2 - dt$mid1)
  
  dt$sample <- as.factor(sample_id)
  dt$resolution <- resolution
  dt$method <- method
  
  dt[, c("sample", "resolution", "method", "chr1", "chr2", "contact_distance_bp")]
}

# Build long table
loop_dists <- data.frame()

for (sample in samples) {
  for (resolution in resolutions) {
    res_kb <- resolution / 1000
    
    FitHiC2_file <- paste0(
      "<path>/FULL_SAMPLES/FITHIC_2026/LOOPS/BEDPES/",
      sample, ".fithic.", resolution, ".loops.bedpe"
    )
    
    hiccups_file <- paste0(
      "<path>/FULL_SAMPLES/loops/",
      "loops_", res_kb, "kb", sample, ".bedpe"
    )
    
    tmp_FitHiC2 <- get_contact_distances(
      file = FitHiC2_file,
      sample_id = sample,
      resolution = resolution,
      method = "FitHiC2",
      header = TRUE
    )
    
    tmp_hiccups <- get_contact_distances(
      file = hiccups_file,
      sample_id = sample,
      resolution = resolution,
      method = "HiCCUPS",
      header = TRUE
    )
    
    loop_dists <- rbind(loop_dists, tmp_FitHiC2, tmp_hiccups)
  }
}

loop_dists$sample <- factor(loop_dists$sample, levels = as.character(samples))

# Plotting function: one plot per method per resolution
make_violin_plot <- function(df, method_name, resolution_value, cap_y = FALSE, ymax = 2e6) {
  p <- ggplot(
    data = df %>% filter(method == method_name, resolution == resolution_value),
    aes(x = sample, y = contact_distance_bp, group = sample)
  ) +
    geom_violin(fill = "lightblue", color = "black", trim = TRUE, scale = "width") +
    geom_boxplot(
      width = 0.12,
      fill = "white",
      color = "black",
      outlier.shape = NA
    ) +
    scale_y_continuous(
      breaks = c(0, 5e5, 1e6, 1.5e6, 2e6),
      labels = c(
        "0",
        expression(5 %*% 10^5),
        expression(1 %*% 10^6),
        expression(1.5 %*% 10^6),
        expression(2 %*% 10^6)
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    xlab("Sample ID") +
    ylab("Contact Distance (bp)") +
    ggtitle(paste0(method_name, " Loop Contact Distance Distributions, ", resolution_value / 1000, " kb"))
  
  if (cap_y) {
    p <- p + coord_cartesian(ylim = c(0, ymax))
  }
  
  p
}

# Generate separate plots per resolution for FitHiC2
for (res in resolutions) {
  p <- make_violin_plot(
    df = loop_dists,
    method_name = "FitHiC2",
    resolution_value = res,
    cap_y = cap_y_at_2e6,
    ymax = y_cap
  )
  
  suffix <- ifelse(cap_y_at_2e6, "_ymax2e6", "")
  ggsave(
    filename = paste0("../QCFIGS/FitHiC2_contact_distance_", res / 1000, "kb_violin", suffix, ".pdf"),
    plot = p,
    width = 12,
    height = 6
  )
}

# Generate separate plots per resolution for HiCCUPS
for (res in resolutions) {
  p <- make_violin_plot(
    df = loop_dists,
    method_name = "HiCCUPS",
    resolution_value = res,
    cap_y = cap_y_at_2e6,
    ymax = y_cap
  )
  
  suffix <- ifelse(cap_y_at_2e6, "_ymax2e6", "")
  ggsave(
    filename = paste0("../QCFIGS/hiccups_contact_distance_", res / 1000, "kb_violin", suffix, ".pdf"),
    plot = p,
    width = 12,
    height = 6
  )
}



##### Get HiCRes results as one panel
library(magick)
library(cowplot)
library(patchwork)

base_dir <- "<path>/FULL_SAMPLES/hicres_results"

dirs <- list.dirs(base_dir, recursive = FALSE)
dirs <- dirs[!grepl("old|QC_|SAMPLE_01|SAMPLE_02|SAMPLE_03|SAMPLE_04|4200", dirs)]

get_sample <- function(path) {
  sub("_HiCRes_outputs$", "", basename(path))
}

plots <- list()

for (d in dirs) {
  img_path <- file.path(d, "prediction.png")
  if (!file.exists(img_path)) next
  
  sample <- get_sample(d)
  
  img <- image_read(img_path)
  img <- image_trim(img)
  
  p <- ggdraw() +
    draw_label(
      sample,
      x = 0.5, y = 0.95,
      hjust = 0.5, vjust = 1,
      fontface = "bold",
      size = 12
    ) +
    draw_image(
      img,
      x = 0, y = 0,
      width = 1, height = 0.98,
      hjust = 0, vjust = 0
    )
  
  plots[[sample]] <- p
}

combined <- wrap_plots(plots, ncol = 4)

ggsave(
  filename = "../QCFIGS/hicres_combined_panels.pdf",
  plot = combined,
  width = 16,
  height = 20
)

##### Attempt to make it in a single figure:
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

all_preds <- data.frame()

for (d in dirs) {
  pred_file <- file.path(d, "prediction.txt")
  if (!file.exists(pred_file)) next
  
  sample <- get_sample(d)
  tmp <- fread(pred_file)
  
  # enforce names explicitly
  colnames(tmp) <- c("Valid_Reads", "Predicted_Resolution")
  
  tmp$sample <- sample
  tmp$Valid_Reads_M <- tmp$Valid_Reads / 1e6
  
  all_preds <- rbind(all_preds, tmp)
}

all_preds$sample <- factor(all_preds$sample, levels = sort(unique(all_preds$sample)))

p_overlay <- ggplot(
  all_preds,
  aes(x = Valid_Reads_M, y = Predicted_Resolution, color = sample, group = sample)
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  theme_bw() +
  xlab("Valid read pairs (M)") +
  ylab("Predicted resolution (bp)") +
  ggtitle("HiCRes Predicted Resolution Curves Across Samples") +
  scale_y_continuous(labels = comma_format()) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) + coord_cartesian(ylim=c(0, 500000))

ggsave(
  "../QCFIGS/hicres_overlay_all_samples_from_prediction_txt.pdf",
  plot = p_overlay,
  width = 12,
  height = 8
)

p_overlay_log <- ggplot(
  all_preds,
  aes(x = Valid_Reads_M, y = Predicted_Resolution, color = sample, group = sample)
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  theme_bw() +
  xlab("Valid read pairs (M)") +
  ylab("Predicted resolution (bp, log scale)") +
  ggtitle("HiCRes Predicted Resolution Curves Across Samples") +
  scale_y_log10(
    breaks = c(2e3, 5e3, 1e4, 2e4, 1e5, 5e5),
    labels = comma_format()
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(1e3, 5e5))
ggsave(
  "../QCFIGS/hicres_overlay_all_samples_from_prediction_txt_logY.pdf",
  plot = p_overlay_log,
  width = 12,
  height = 8
)


##### Take a look at unique valid pairs and inferred resolution:
qc_base_dir <- file.path(base_dir, "QC_outputs")
# -------------------------
# Helper functions
# -------------------------

read_valid_pairs <- function(sample_id) {
  qc_dir <- file.path(qc_base_dir, paste0(sample_id, "_QC_outputs"))
  qc_file <- list.files(qc_dir, pattern = "^hicup_filter_summary_.*\\.txt$", full.names = TRUE)
  
  if (length(qc_file) != 1) {
    stop(paste("Expected exactly 1 HiCUP summary file for sample", sample_id,
               "but found", length(qc_file)))
  }
  
  dt <- fread(qc_file)
  if (!"Valid_pairs" %in% colnames(dt)) {
    stop(paste("Valid_pairs column not found in", qc_file))
  }
  
  data.frame(
    sample = as.character(sample_id),
    Valid_pairs = as.numeric(dt$Valid_pairs[1])
  )
}

read_prediction <- function(sample_id) {
  pred_file <- file.path(
    base_dir,
    paste0(sample_id, "_HiCRes_outputs"),
    "prediction.txt"
  )
  
  if (!file.exists(pred_file)) {
    stop(paste("Missing prediction.txt for sample", sample_id))
  }
  
  dt <- fread(pred_file)
  colnames(dt) <- c("Valid_Reads", "Predicted_Resolution")
  dt$sample <- as.character(sample_id)
  dt
}

# nearest predicted resolution for actual valid-pair count
get_resolution_at_sample_depth <- function(pred_dt, valid_pairs) {
  pred_dt <- pred_dt %>%
    mutate(diff = abs(Valid_Reads - valid_pairs)) %>%
    arrange(diff)
  
  pred_dt$Predicted_Resolution[1]
}

# optional linear interpolation if actual valid_pairs falls between grid points
get_resolution_interpolated <- function(pred_dt, valid_pairs) {
  approx(
    x = pred_dt$Valid_Reads,
    y = pred_dt$Predicted_Resolution,
    xout = valid_pairs,
    rule = 2
  )$y
}

# -------------------------
# Build summary tables
# -------------------------

valid_pairs_df <- do.call(
  rbind,
  lapply(samples, read_valid_pairs)
)

prediction_list <- lapply(samples, read_prediction)
all_preds <- do.call(rbind, prediction_list)

resolution_df <- lapply(samples, function(s) {
  s_char <- as.character(s)
  
  vp <- valid_pairs_df %>%
    filter(sample == s_char) %>%
    pull(Valid_pairs)
  
  pred_dt <- all_preds %>%
    filter(sample == s_char)
  
  # choose one:
  # predicted_res <- get_resolution_at_sample_depth(pred_dt, vp)     # nearest point
  predicted_res <- get_resolution_interpolated(pred_dt, vp)          # smoother / better
  
  data.frame(
    sample = s_char,
    Valid_pairs = vp,
    Predicted_Resolution = predicted_res
  )
}) %>% bind_rows()

# preserve sample order
valid_pairs_df$sample <- factor(valid_pairs_df$sample, levels = as.character(samples))
resolution_df$sample <- factor(resolution_df$sample, levels = as.character(samples))

# -------------------------
# Plot 1: valid pairs barplot
# -------------------------

p_valid_pairs <- ggplot(valid_pairs_df, aes(x = sample, y = Valid_pairs)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  xlab("Sample ID") +
  ylab("Valid pairs") +
  ggtitle("Valid Pairs per Sample") +
  scale_y_continuous(labels = comma_format(), breaks=c(250000000, 500000000, 750000000, 1000000000, 1250000000))

ggsave(
  "../QCFIGS/hicup_valid_pairs_per_sample.pdf",
  plot = p_valid_pairs,
  width = 10,
  height = 6
)

# -------------------------
# Plot 2: inferred HiCRes resolution barplot
# -------------------------

p_resolution <- ggplot(resolution_df, aes(x = sample, y = Predicted_Resolution)) +
  geom_col(fill = "firebrick2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  xlab("Sample ID") +
  ylab("Predicted resolution (bp)") +
  ggtitle("HiCRes Predicted Resolution at Observed Valid-Pair Depth") +
  scale_y_continuous(labels = comma_format(), breaks=c(1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3)) + coord_cartesian(ylim=c(0, 10001))

ggsave(
  "../QCFIGS/hicres_predicted_resolution_at_observed_depth.pdf",
  plot = p_resolution,
  width = 10,
  height = 6
)

# -------------------------
# Optional: merged summary table
# -------------------------

summary_df <- valid_pairs_df %>%
  left_join(resolution_df %>% select(sample, Predicted_Resolution), by = "sample")

write.table(
  summary_df,
  file = "../QCFIGS/hicres_valid_pairs_and_resolution_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


###### Sequence duplication levels from multiQC
base_dir <- "<path>/FULL_SAMPLES"

# Sample order to use in the plot
sample_order <- samples

# Helper: read one sample's multiqc duplication table
read_duplication_file <- function(sample_id, base_dir) {
  infile <- file.path(
    base_dir,
    as.character(sample_id),
    "multiqc",
    "mqc_fastqc_sequence_counts_plot_1.txt"
  )
  
  if (!file.exists(infile)) {
    warning(paste("Missing file for sample", sample_id, ":", infile))
    return(NULL)
  }
  
  dt <- fread(infile)
  
  required_cols <- c("Sample", "Unique Reads", "Duplicate Reads")
  if (!all(required_cols %in% colnames(dt))) {
    stop(
      paste(
        "File", infile, "is missing required columns. Found:",
        paste(colnames(dt), collapse = ", ")
      )
    )
  }
  
  # Sum across read1/read2 rows
  unique_reads <- sum(as.numeric(dt[["Unique Reads"]]), na.rm = TRUE)
  duplicate_reads <- sum(as.numeric(dt[["Duplicate Reads"]]), na.rm = TRUE)
  total_reads <- unique_reads + duplicate_reads
  
  pct_duplicate <- if (total_reads > 0) {
    100 * duplicate_reads / total_reads
  } else {
    NA_real_
  }
  
  data.frame(
    sample = as.character(sample_id),
    unique_reads = unique_reads,
    duplicate_reads = duplicate_reads,
    total_reads = total_reads,
    pct_duplicate = pct_duplicate,
    stringsAsFactors = FALSE
  )
}

# Build summary table
dup_df <- bind_rows(lapply(sample_order, read_duplication_file, base_dir = base_dir))

# Preserve plotting order
dup_df$sample <- factor(dup_df$sample, levels = as.character(sample_order))

# Write summary table
write.table(
  dup_df,
  file = "../QCFIGS/percent_duplicate_reads_per_sample.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Plot: percent duplicate reads per sample
p_dup <- ggplot(dup_df, aes(x = sample, y = pct_duplicate)) +
  geom_col(fill = "steelblue") +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_duplicate)),
    vjust = -0.25,
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  xlab("Sample ID") +
  ylab("Duplicate reads (%)") +
  ggtitle("Percent Duplicate Reads per Sample (Raw FastQ)") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(ylim = c(0, max(dup_df$pct_duplicate, na.rm = TRUE) * 1.08))

ggsave(
  "../QCFIGS/percent_duplicate_reads_per_sample.pdf",
  plot = p_dup,
  width = 11,
  height = 6
)


###### HiC-Pro Contact Statistics BarPlots
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

base_dir <- "<path>/FULL_SAMPLES"

sample_order <- samples

contact_categories <- c(
  "Unique: cis <= 20Kbp",
  "Unique: cis > 20Kbp",
  "Unique: trans",
  "Duplicate read pairs"
)

contact_colors <- c(
  "Unique: cis <= 20Kbp" = "#143fd6",  # dark blue
  "Unique: cis > 20Kbp" = "#7893e2",   # light blue
  "Unique: trans" = "#06a02d",         # green
  "Duplicate read pairs" = "#ada7a7"   # gray
)

read_contact_file <- function(sample_id, base_dir) {
  infile <- file.path(
    base_dir,
    as.character(sample_id),
    "multiqc",
    "mqc_hicpro_contact_plot_1.txt"
  )
  
  if (!file.exists(infile)) {
    warning(paste("Missing file for sample", sample_id, ":", infile))
    return(NULL)
  }
  
  dt <- fread(infile)
  
  required_cols <- c(
    "Sample",
    "Unique: cis <= 20Kbp",
    "Unique: cis > 20Kbp",
    "Unique: trans",
    "Duplicate read pairs"
  )
  
  if (!all(required_cols %in% colnames(dt))) {
    stop(
      paste(
        "File", infile, "is missing required columns. Found:",
        paste(colnames(dt), collapse = ", ")
      )
    )
  }
  
  # Expect one row per sample file; if multiple rows exist, sum them
  out <- data.frame(
    sample = as.character(sample_id),
    `Unique: cis <= 20Kbp` = sum(as.numeric(dt[["Unique: cis <= 20Kbp"]]), na.rm = TRUE),
    `Unique: cis > 20Kbp`  = sum(as.numeric(dt[["Unique: cis > 20Kbp"]]), na.rm = TRUE),
    `Unique: trans`        = sum(as.numeric(dt[["Unique: trans"]]), na.rm = TRUE),
    `Duplicate read pairs` = sum(as.numeric(dt[["Duplicate read pairs"]]), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out
}

# Read all sample files
contact_df_wide <- bind_rows(
  lapply(sample_order, read_contact_file, base_dir = base_dir)
)

# Preserve order
contact_df_wide$sample <- factor(contact_df_wide$sample, levels = as.character(sample_order))

# Fix colnames
colnames(contact_df_wide)[2:5] <- c("Unique: cis <= 20Kbp", "Unique: cis > 20Kbp", "Unique: trans", "Duplicate read pairs")

# Write summary table
write.table(
  contact_df_wide,
  file = "../QCFIGS/hicpro_contact_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Convert to long format
contact_df_long <- contact_df_wide %>%
  pivot_longer(
    cols = all_of(contact_categories),
    names_to = "category",
    values_to = "count"
  )

contact_df_long$category <- factor(contact_df_long$category, levels = contact_categories)

# Add percentages
contact_df_long <- contact_df_long %>%
  group_by(sample) %>%
  mutate(
    total = sum(count, na.rm = TRUE),
    percent = 100 * count / total
  ) %>%
  ungroup()

contact_df_long$category <- factor(
  contact_df_long$category,
  levels = rev(contact_categories)
)

# -------------------------
# Plot 1: stacked counts
# -------------------------
p_counts <- ggplot(contact_df_long, aes(x = sample, y = count, fill = category)) +
  geom_col() +
  scale_fill_manual(values = contact_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  xlab("Sample ID") +
  ylab("Read pair count") +
  ggtitle("HiC-Pro Contact Statistics per Sample") +
  scale_y_continuous(labels = comma_format()) +
  guides(fill = guide_legend(title = NULL))

ggsave(
  "../QCFIGS/hicpro_contact_statistics_counts_stacked.pdf",
  plot = p_counts,
  width = 12,
  height = 6
)

# -------------------------
# Plot 2: stacked percentages
# -------------------------
p_percent <- ggplot(contact_df_long, aes(x = sample, y = percent, fill = category)) +
  geom_col() +
  scale_fill_manual(values = contact_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  xlab("Sample ID") +
  ylab("Percentage of read pairs") +
  ggtitle("HiC-Pro Contact Statistics per Sample (%)") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(fill = guide_legend(title = NULL))

ggsave(
  "../QCFIGS/hicpro_contact_statistics_percent_stacked.pdf",
  plot = p_percent,
  width = 12,
  height = 6
)


##### Patchwork together the coverage and saddle plots, note usage of convert from command line for this (i.e. the below are bash scripts, NOT R code!):
### COVERAGE PLOTS
Decided to re-make these from scratch in python to ensure they are connected and file is not HUGE


### SADDLE PLOTS (done in <compute_environment>, so I could install pdfjam)
cd <path>/FULL_SAMPLES/QCFIGS

pdfjam [0-9]*.1000000_saddle_plot.pdf \
--nup 5x3 \
--landscape \
--paper a1paper \
--outfile all_samples_saddle_plot_panel.pdf
