# visualize_scores.R
# -------------------
# Snakemake-integrated visualisation of splice junction scoring results
# produced by score_splice_junctions.py.
#
# Plots produced (one page each in the output PDF):
#   1. Per-junction heatmap   — match_status per tool, faceted by event_type GROUP
#   2. Junction count bar     — exact/partial/fp stacked per tool, n_truth line
#   2b. Recall vs Precision scatter
#   2C-i.  Stacked bar: TP per tool, coloured by event type
#   2C-ii. Grouped bar: TP per event type, bars per tool
#   A1. Arrow scatter: individual → individual+stringtie
#   A2. Arrow scatter: each individual tool → consensus
#   B.  Dumbbell plot: individual vs +stringtie per pair
#   C.  Delta bar chart: Δrecall and Δprecision per pair
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggrepel)   # repel overlapping text labels apart automatically
})

# ---------------------------------------------------------------------------
# 0. Load inputs
# ---------------------------------------------------------------------------
cat(sprintf("Loading data...\n"))
summary_df <- read_tsv(snakemake@input[["summary"]],           show_col_types = FALSE)
junc_ov    <- read_tsv(snakemake@input[["junction_overview"]], show_col_types = FALSE)
log_file    <- snakemake@log[[1]]
out_pdf    <- snakemake@output[["plots"]]

tool_order <- snakemake@params[["tool_names"]]
strict <- snakemake@params[["strict"]]

summary_df <- summary_df %>% mutate(tool = factor(tool, levels = tool_order))
cat(sprintf("Initiating logging...\n"))
cat("Preparing Visualization...")

fontfamily <- "Times"
fontsize_title <- 30
base_size <- 25
p4 <- NULL

majiq_color <- "#cc79a7"
whippet_color <- "#d55e00"
leafcutter_color <- "#0072b2"
spliceat_color <- "#000000"

# ---------------------------------------------------------------------------
# 1. Event-type grouping
#
# Patterns are case-insensitive and applied in order; the first match wins.
# Anything not matched by the first three rules keeps its original name as
# its own group. Adjust the regex strings here if your naming differs.
# ---------------------------------------------------------------------------

assign_group <- function(event_type) {
  case_when(
    str_detect(event_type, regex("alternative_donor_early",   ignore_case = TRUE)) ~ "cADus",
    str_detect(event_type, regex("alternative_donor_late",   ignore_case = TRUE)) ~ "cADds",
    str_detect(event_type, regex("alternative_acceptor_early",ignore_case = TRUE)) ~ "cAAus",
    str_detect(event_type, regex("alternative_acceptor_late",ignore_case = TRUE)) ~ "cAAds",
    str_detect(event_type, regex("ES|exon_skip",        ignore_case = TRUE)) ~ "cES",
    str_detect(event_type, regex("IR|intron_retention",  ignore_case = TRUE)) ~ "cIR",
    TRUE ~ event_type   # all other event types keep their name as group
  )
}

# Fixed display order for known groups; any novel groups are appended after
known_group_order <- c(
  "cADus",
  "cADds",
  "cAAus",
  "cAAds",
  "cES",
  "cIR"
)


# ---------------------------------------------------------------------------

FS_BASE    <- 28          # base_size for theme_bw  → axes, ticks, strip text
FS_LABEL   <- 32          # axis titles
FS_TICK    <- 26          # tick labels
FS_LEGEND  <- 26          # legend title + keys
FS_REPEL   <- 7           # geom_text_repel label size (in mm-ish ggplot units)


# Build the full group-level order (known groups first, then any extras)
all_groups <- summary_df %>%
  filter(event_type != "ALL") %>%
  distinct(event_type) %>%
  mutate(group = assign_group(event_type)) %>%
  pull(group) %>%
  unique()

group_order <- c(
  known_group_order[known_group_order %in% all_groups],
  setdiff(all_groups, known_group_order)
)

# ---------------------------------------------------------------------------
# Shared colour palette
# ---------------------------------------------------------------------------

status_colors <- c(
  exact   = "#49714a",
  partial = "#FFC107",
  fn      = "#af0000",
  fp      = "#ff0000"
)

# ---------------------------------------------------------------------------
# Plot 1 — Per-junction heatmap
# Source : junction_overview.tsv
#
# Rows  = individual truth junctions (y-labels hidden)
# Cols  = tools
# Fill  = match_status (exact / partial / fn)
# Facet = event_type GROUP (dynamically grouped)
# ---------------------------------------------------------------------------

meta_cols <- c("event_id", "event_type")
tool_cols  <- setdiff(names(junc_ov), meta_cols)

junc_ov <- junc_ov %>%
  mutate(group = assign_group(event_type)) %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group, event_type, event_id) %>%
  mutate(junction_id = row_number())

junc_long <- junc_ov %>%
  pivot_longer(
    cols      = all_of(tool_cols),
    names_to  = "tool",
    values_to = "match_status"
  ) %>%
  mutate(
    tool         = factor(tool, levels = tool_order),
    match_status = factor(match_status, levels = c("exact", "partial", "fn")),
    junction_id  = fct_rev(factor(junction_id))
  )
cat("Plotting...")

p1 <- ggplot(junc_long, aes(x = tool, y = junction_id, fill = match_status)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values   = status_colors,
    na.value = "#E0E0E0",
    name     = "Match status",
    drop     = FALSE
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  labs(x = "Tool", y = NULL) +
  theme_minimal(base_size = base_size) +
  theme(
    text             = element_text(family = fontfamily),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = FS_TICK),
    axis.title.x     = element_text(size = FS_LABEL),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0, "lines"),   # ← removes inter-facet gap
    strip.text.y     = element_text(angle = 0, hjust = 0, face = "bold",
                                    size = FS_LABEL),
    strip.clip       = "off",              # ← prevents strip labels being clipped
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.title     = element_text(size = FS_LEGEND, face = "bold"),
    legend.text      = element_text(size = FS_LEGEND),
    legend.key.size  = unit(1.4, "lines"),
    plot.title       = element_blank()
  )

# ---------------------------------------------------------------------------
# Plot 2 — Stacked bar: predicted junctions as exact / partial / fp per tool
# Source : summary.tsv (event_type == "ALL")
# ---------------------------------------------------------------------------

sum_all <- summary_df %>% filter(event_type == "ALL")

bar2_long <- sum_all %>%
  select(tool, exact, partial, fp) %>%
  pivot_longer(
    cols      = c(exact, partial, fp),
    names_to  = "category",
    values_to = "count"
  ) %>%
  mutate(category = factor(category, levels = c("fp", "partial", "exact"))) # exact at bottom

n_truth_val <- sum_all %>% slice(1) %>% pull(n_truth)

p2 <- ggplot(bar2_long, aes(x = tool, y = count, fill = category)) +
  geom_col(position = "stack", width = 0.6) +
  geom_hline(
    yintercept = n_truth_val,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.7
  ) +
  annotate(
    "text",
    x     = length(tool_order) + 0.45,
    y     = n_truth_val,
    label = paste0("n_truth = ", n_truth_val),
    hjust = 1, vjust = -0.45,
    size  = 5
  ) +
  scale_fill_manual(values = status_colors, name = "Category") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Predicted Junctions by Match Category per Tool",
    x     = "Tool",
    y     = "Number of predicted junctions"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    text            = element_text(family = fontfamily),
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = fontsize_title)
  )

if (strict) {
  recall_precision <- summary_df %>%
  filter(event_type == "ALL") %>%
  select(tool, exact, partial, fp, n_truth) %>%
  mutate(
    recall    = (exact) / n_truth,
    precision = (exact) / (exact + partial + fp)
  ) %>%
  select(tool, recall, precision)
} else {
  recall_precision <- summary_df %>%
    filter(event_type == "ALL") %>%
    select(tool, exact, partial, fp, n_truth) %>%
    mutate(
      recall    = (exact + partial) / n_truth,
      precision = (exact + partial) / (exact + partial + fp)
  ) %>%
  select(tool, recall, precision)
}

# ---------------------------------------------------------------------------
# Plot 2b — Scatter plot of recall vs precision
# ---------------------------------------------------------------------------

p2b <- ggplot(recall_precision, aes(x = recall, y = precision, label = tool)) +
  geom_point(size = 3) +
  geom_text_repel(
    size          = 4,
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = Inf,
    direction     = "both",
    seed          = 42
  ) +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(NA, 1), ylim = c(NA, 1), clip = "off") +
  labs(
    title = "Recall vs Precision by Tool and Combination",
    x     = "Recall",
    y     = "Precision"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    aspect.ratio    = 1,
    text            = element_text(family = fontfamily),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = fontsize_title)
  )

# ---------------------------------------------------------------------------
# Plot 2C — True positives broken down by event type
#
# Two complementary views of the same data (exact + partial per tool per group):
#   2C-i  stacked bar  — one bar per tool, segments = event type groups
#                        → easy to compare total TP height across tools
#   2C-ii grouped bar  — one cluster per event type, bars = tools
#                        → easy to compare tools within a single event type
# ---------------------------------------------------------------------------

tp_by_type <- summary_df %>%
  filter(event_type != "ALL") %>%
  mutate(
    group = assign_group(event_type),
    # Reverse group order for 2C-i so the first group sits at the top of the stack
    group = factor(group, levels = rev(group_order)),
    tool  = factor(tool,  levels = tool_order),
    tp    = if (strict) exact else exact + partial
  )

tp_label <- if (strict) "True positives (exact only)" else "True positives (exact + partial)"

# Colour palette for event type groups — keyed on rev(group_order) so colours
# are consistent with the legend even after reversing the stack order
group_colors <- setNames(
  hue_pal()(length(group_order)),
  rev(group_order)
)

# Colour palette for tools (one hue per tool, reused in 2C-ii)
tool_colors <- setNames(
  hue_pal()(length(tool_order)),
  tool_order
)

# Fixed y-axis: breaks every 30, hard ceiling at 180
y_scale_tp <- scale_y_continuous(
  breaks = seq(0, 180, by = 30),
  limits = c(0, 180),
  expand = expansion(mult = c(0, 0))
)
y_scale_tp_2 <- scale_y_continuous(
  breaks = seq(0, 30, by = 5),
  limits = c(0, 30),
  expand = expansion(mult = c(0, 0))
)

# -- 2C-i: stacked bar, one bar per tool ------------------------------------
p2c_stacked <- ggplot(tp_by_type, aes(x = tool, y = tp, fill = group)) +
  geom_col(position = "stack", width = 0.65, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = group_colors, name = "Event type") +
  y_scale_tp +
  labs(
    title = "True Positives per Event Type and Tool",
    x     = "Tool",
    y     = tp_label
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    text               = element_text(family = fontfamily),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = fontsize_title)
  )

# -- 2C-ii: grouped bar, one cluster per event type -------------------------
p2c_grouped <- ggplot(tp_by_type, aes(x = group, y = tp, fill = tool)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tool_colors, name = "Tool") +
  y_scale_tp_2 +
  labs(
    title = "True Positives per Event Type and Tool",
    x     = "Event type",
    y     = tp_label
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    text               = element_text(family = fontfamily),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = fontsize_title)
  )

# ---------------------------------------------------------------------------
# Shared pairing structure
#
# tool_order is assumed to contain an even number of entries organised as
# consecutive pairs:  (individual, individual+stringtie, ...)
# The last pair is the consensus tool.
# ---------------------------------------------------------------------------

n_tools  <- length(tool_order)
n_pairs  <- n_tools / 2

pairs_df <- tibble(
  pair_idx     = seq_len(n_pairs),
  individual   = tool_order[seq(1, n_tools, by = 2)],
  stringtie    = tool_order[seq(2, n_tools, by = 2)],
  is_consensus = c(rep(FALSE, n_pairs - 1), TRUE),
  pair_label   = tool_order[seq(1, n_tools, by = 2)]
)

# Attach recall/precision coordinates for both members of each pair
arrow_df <- pairs_df %>%
  left_join(recall_precision, by = c("individual" = "tool")) %>%
  rename(recall_ind = recall, precision_ind = precision) %>%
  left_join(recall_precision, by = c("stringtie" = "tool")) %>%
  rename(recall_st = recall, precision_st = precision)

# Consensus coordinates (single reference point used in plot A2)
consensus_row  <- pairs_df %>% filter(is_consensus)
consensus_ind  <- recall_precision %>% filter(tool == consensus_row$individual)
consensus_st   <- recall_precision %>% filter(tool == consensus_row$stringtie)

# Non-consensus pairs only (used in plots A2, B, C)
pairs_noncons  <- pairs_df  %>% filter(!is_consensus)
arrow_noncons  <- arrow_df  %>% filter(!is_consensus)

# Shared arrow styling
arr <- arrow(length = unit(0.25, "cm"), type = "closed")

# Colour: one hue per non-consensus pair, consensus always black
pair_colors <- setNames(
  c(leafcutter_color, majiq_color, whippet_color, spliceat_color),
  pairs_df$pair_label
)

# ---------------------------------------------------------------------------
# Plot A1 — Arrow scatter: individual → individual + stringtie
# legend: top-left INSIDE the coordinate system
# ---------------------------------------------------------------------------

label_a1 <- bind_rows(
  arrow_df %>% transmute(
    label      = individual,
    recall     = recall_ind,
    precision  = precision_ind,
    pair_label
  ),
  arrow_df %>% transmute(
    label      = stringtie,
    recall     = recall_st,
    precision  = precision_st,
    pair_label
  )
)

pA1 <- ggplot() +
  geom_segment(
    data = arrow_df,
    aes(x     = recall_ind, y    = precision_ind,
        xend  = recall_st,  yend = precision_st,
        color = pair_label),
    arrow     = arr,
    linewidth = 0.9,
    lineend   = "round"
  ) +
  geom_point(
    data  = label_a1,
    aes(x = recall, y = precision, color = pair_label),
    size  = 3
  ) +
  geom_text_repel(
    data          = label_a1,
    aes(x = recall, y = precision, label = label, color = pair_label),
    size          = FS_REPEL,
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = Inf,
    direction     = "both",
    seed          = 42,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = pair_colors, name = "Pair") +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(NA, 1), ylim = c(NA, 1)) +
  labs(x = "Recall", y = "Precision") +
  theme_bw(base_size = FS_BASE) +
  theme(
    aspect.ratio    = 1,
    text            = element_text(family = fontfamily),
    plot.title      = element_blank(),

    # axis labels
    axis.title      = element_text(size = FS_LABEL),
    axis.text       = element_text(size = FS_TICK),

    # legend INSIDE — top-left corner
    legend.position        = c(0.98, 0.02),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 0.7),
                                          colour = "grey70", linewidth = 0.4),
    legend.title           = element_text(size = FS_LEGEND, face = "bold"),
    legend.text            = element_text(size = FS_LEGEND),
    legend.key.size        = unit(1.4, "lines"),
    legend.margin          = margin(4, 6, 4, 6),
  )

# ---------------------------------------------------------------------------
# Plot A2 — Arrow scatter: each individual tool → consensus
# legend: top (outside, horizontal) — consensus point has no colour so inside
#         would be cluttered; top is cleaner
# ---------------------------------------------------------------------------

a2_ind <- arrow_noncons %>%
  transmute(
    label      = individual,
    x          = recall_ind,
    y          = precision_ind,
    xend       = consensus_ind$recall,
    yend       = consensus_ind$precision,
    pair_label,
    line_type  = "solid"
  )

a2_st <- arrow_noncons %>%
  transmute(
    label      = stringtie,
    x          = recall_st,
    y          = precision_st,
    xend       = consensus_st$recall,
    yend       = consensus_st$precision,
    pair_label,
    line_type  = "dashed"
  )

a2_arrows <- a2_ind

a2_points <- bind_rows(
  a2_arrows %>% transmute(label, recall = x, precision = y, pair_label),
  tibble(
    label      = c(consensus_ind$tool, consensus_st$tool),
    recall     = c(consensus_ind$recall, consensus_st$recall),
    precision  = c(consensus_ind$precision, consensus_st$precision),
    pair_label = "Consensus"
  )
)

pA2 <- ggplot() +
  geom_segment(
    data = a2_arrows,
    aes(x = x, y = y, xend = xend, yend = yend,
        color = pair_label, linetype = line_type),
    arrow     = arr,
    linewidth = 0.9,
    lineend   = "round"
  ) +
  geom_point(
    data  = a2_points %>% filter(pair_label == "Consensus"),
    aes(x = recall, y = precision),
    shape = 18, size = 6, color = "black"
  ) +
  geom_point(
    data  = a2_points %>% filter(pair_label != "Consensus"),
    aes(x = recall, y = precision, color = pair_label),
    size  = 3
  ) +
  geom_text_repel(
    data          = a2_points,
    aes(x = recall, y = precision, label = label, color = pair_label),
    size          = FS_REPEL,
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = Inf,
    direction     = "both",
    seed          = 42,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = pair_colors, name = "Tool") +
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(NA, 1), ylim = c(NA, 1)) +
  labs(x = "Recall", y = "Precision") +
  theme_bw(base_size = FS_BASE) +
  theme(
    aspect.ratio    = 1,
    text            = element_text(family = fontfamily),
    plot.title      = element_blank(),

    # axis labels
    axis.title      = element_text(size = FS_LABEL),
    axis.text       = element_text(size = FS_TICK),

    # legend at top, horizontal
    legend.position        = c(0.98, 0.02),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 0.7),
                                          colour = "grey70", linewidth = 0.4),
    legend.title           = element_text(size = FS_LEGEND, face = "bold"),
    legend.text            = element_text(size = FS_LEGEND),
    legend.key.size        = unit(1.4, "lines"),
    legend.margin          = margin(4, 6, 4, 6),
    legend.box.margin      = margin(0, 0, -6, 0),   # pull legend closer to plot
  )
# ---------------------------------------------------------------------------
# Plot B — Dumbbell plot: individual vs +stringtie per pair
# ---------------------------------------------------------------------------

dumbbell_df <- arrow_df %>%
  transmute(
    pair_label,
    is_consensus,
    Recall    = recall_ind,    Recall_st    = recall_st,
    Precision = precision_ind, Precision_st = precision_st
  ) %>%
  pivot_longer(
    cols            = c(Recall, Precision, Recall_st, Precision_st),
    names_to        = "metric_raw",
    values_to       = "value"
  ) %>%
  mutate(
    metric  = if_else(str_detect(metric_raw, "Recall"), "Recall", "Precision"),
    variant = if_else(str_detect(metric_raw, "_st"),    "Individual + StringTie", "Individual")
  ) %>%
  select(-metric_raw) %>%
  mutate(pair_label = factor(pair_label, levels = rev(pairs_df$pair_label)))

pB <- ggplot(dumbbell_df, aes(y = pair_label)) +
  geom_line(
    aes(x = value, group = interaction(pair_label, metric),
        color = pair_label),
    linewidth = 1.2, alpha = 0.5
  ) +
  geom_point(
    aes(x = value, shape = variant, color = pair_label),
    size = 4
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype   = "dashed",
    color      = "grey60",
    linewidth  = 0.6
  ) +
  scale_color_manual(values = pair_colors, guide = "none") +
  scale_shape_manual(
    values = c("Individual" = 19, "Individual + StringTie" = 17),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(add = c(0.02, 0.02))) +
  facet_wrap(~ metric, ncol = 2) +
  labs(
    title = "Effect of StringTie Filtering on Recall and Precision",
    x     = "Value",
    y     = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    text            = element_text(family = fontfamily),
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = fontsize_title),
    plot.subtitle   = element_text(size = base_size - 4),
    strip.text      = element_text(face = "bold")
  )

# ---------------------------------------------------------------------------
# Plot C — Delta bar chart: Δrecall and Δprecision per pair
# ---------------------------------------------------------------------------

delta_df <- arrow_df %>%
  transmute(
    pair_label,
    is_consensus,
    delta_Recall    = recall_st    - recall_ind,
    delta_Precision = precision_st - precision_ind
  ) %>%
  pivot_longer(
    cols      = c(delta_Recall, delta_Precision),
    names_to  = "metric",
    values_to = "delta"
  ) %>%
  mutate(
    metric     = str_remove(metric, "delta_"),
    pair_label = factor(pair_label, levels = pairs_df$pair_label),
    direction  = if_else(delta >= 0, "Increase", "Decrease")
  )

pC <- ggplot(delta_df, aes(x = pair_label, y = delta, fill = direction)) +
  geom_col(
    aes(alpha = is_consensus),
    width = 0.6, color = "grey30", linewidth = 0.3
  ) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30") +
  geom_vline(
    xintercept = n_pairs - 0.5,
    linetype   = "dashed",
    color      = "grey60",
    linewidth  = 0.6
  ) +
  scale_fill_manual(
    values = c("Increase" = "#4CAF50", "Decrease" = "#F44336"),
    name   = NULL
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.65),
    guide  = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  facet_wrap(~ metric, ncol = 2) +
  labs(
    title = "Change in Recall and Precision upon StringTie Filtering",
    x     = NULL,
    y     = "Change in Percentage Points"
  ) +
  theme_bw(base_size = 18) +
  theme(
    text            = element_text(family = fontfamily),
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(face = "bold", size = fontsize_title),
    plot.subtitle   = element_text(size = base_size - 4),
    strip.text      = element_text(face = "bold")
  )

# ---------------------------------------------------------------------------
# Write all plots to a multi-page PDF
# ---------------------------------------------------------------------------

n_junctions <- nrow(junc_ov)
p1_height   <- max(7, min(24, 3 + n_junctions * 0.12))

pdf(out_pdf, width = 10, height = 10)
print(p1 + theme(plot.margin = margin(10, 10, 10, 10)))
print(p2)
print(p2b)
print(p2c_stacked)
print(p2c_grouped)
print(pA1)
print(pA2)
print(pB)
print(pC)
dev.off()

message("Plots written to: ", out_pdf)