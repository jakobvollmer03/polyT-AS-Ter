suppressMessages({
    library(data.table)
    library(dplyr)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(biomaRt)
    library(ggplot2)
   # library(showtext)
})


input_dir  <- "/mnt/gtklab01/dbsv771/iclip/prep"
output_dir <- "/mnt/gtklab01/dbsv771/iclip"
clip_table <- "/mnt/gtklab01/dbsv771/iclip/sumgene.tab"
background_gtf <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r0_get_ref/Mus_musculus/Mus_musculus_GRCm39_115_chr_filtered.gtf"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------
# Tool name helper
# More-specific patterns (combined tools) must come before
# their component patterns to avoid premature matching.
# --------------------------------------------------

rename_tool <- function(raw_name) {
    case_when(
        grepl("whippet_stringtie",   raw_name, ignore.case = TRUE) ~ "W + S",
        grepl("majiq_stringtie",     raw_name, ignore.case = TRUE) ~ "M + S",
        grepl("leafcutter_stringtie",raw_name, ignore.case = TRUE) ~ "L + S",
        grepl("ds_consensus",        raw_name, ignore.case = TRUE) ~ "Cons",
        grepl("spliceat",            raw_name, ignore.case = TRUE) ~ "SpliCeAT",
        grepl("whippet",             raw_name, ignore.case = TRUE) ~ "W",
        grepl("majiq",               raw_name, ignore.case = TRUE) ~ "M",
        grepl("leafcutter",          raw_name, ignore.case = TRUE) ~ "L",
        TRUE ~ raw_name
    )
}

# --------------------------------------------------
# Load iCLIP once
# --------------------------------------------------

clip_df <- read.table(
    clip_table,
    skip = 1,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
) %>%
    filter(
        gene_segments_biotypes == "protein_coding",
        total_same_count_sum >= 1
    )

iclip_symbols <- unique(clip_df$gene_name)

map_to_ensembl <- function(genes) {

    map1 <- suppressWarnings(
        bitr(
            genes,
            fromType = "SYMBOL",
            toType = "ENSEMBL",
            OrgDb = org.Mm.eg.db
        )
    )

    map2 <- suppressWarnings(
        bitr(
            genes,
            fromType = "ALIAS",
            toType = "ENSEMBL",
            OrgDb = org.Mm.eg.db
        )
    )

    bind_rows(map1, map2) %>%
        distinct(ENSEMBL, .keep_all = TRUE)
}

iclip_map <- map_to_ensembl(iclip_symbols)

iclip_genes <- unique(iclip_map$ENSEMBL)

cat(sprintf(
    "Mapped %d iCLIP genes\n",
    length(iclip_genes)
))

# --------------------------------------------------
# Background
# --------------------------------------------------

gtf <- fread(
    background_gtf,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    data.table = FALSE
)
background_genes <- gtf %>%
#    filter( # uncomment if gtf is unfiltered
#        V3 == "gene",
#        grepl("protein_coding", V9)
#    ) %>%
    mutate(
        gene_id = sub(".*gene_id \"([^\"]+)\".*", "\\1", V9)
    ) %>%
    pull(gene_id) %>%
    unique()

cat(sprintf(
    "Background gene set contains %d genes\n",
    length(background_genes)
))

# --------------------------------------------------
# Tool files
# --------------------------------------------------

tool_files <- list.files(
    input_dir,
    pattern = "_go_ids\\.txt$",
    full.names = TRUE
)

# --------------------------------------------------
# Results containers
# --------------------------------------------------

summary_results <- list()
gene_level_results <- list()

# --------------------------------------------------
# Iterate tools
# --------------------------------------------------

for (f in tool_files) {

    tool_name <- rename_tool(
        sub("_go_ids\\.txt$", "", basename(f))
    )

    cat(sprintf(
        "\nProcessing %s\n",
        tool_name
    ))

    tool_genes <- fread(f)$gene_id %>%
        unique()

    overlap <- intersect(
        tool_genes,
        iclip_genes
    )

    a <- length(overlap)
    b <- length(tool_genes) - a
    c <- length(iclip_genes) - a
    d <- length(background_genes) - (a + b + c)

    fisher <- fisher.test(
        matrix(
            c(a,b,c,d),
            nrow = 2
        )
    )

    # ---------------------------------------
    # Summary statistics
    # ---------------------------------------

    summary_results[[tool_name]] <- data.frame(
        tool = tool_name,
        n_tool_genes = length(tool_genes),
        n_iclip_genes = length(iclip_genes),
        overlap = a,
        overlap_fraction = a / length(tool_genes),
        odds_ratio = fisher$estimate,
        odds_ratio_low = fisher$conf.int[1],
        odds_ratio_high = fisher$conf.int[2],
        fisher_p = fisher$p.value
    )

    # ---------------------------------------
    # Gene-level table
    # ---------------------------------------

    gene_level_results[[tool_name]] <-
        data.frame(
            gene_id = tool_genes,
            tool = tool_name,
            tdp43_target =
                as.integer(
                    tool_genes %in% iclip_genes
                )
        )
}

# --------------------------------------------------
# Combine results
# --------------------------------------------------

summary_df <- bind_rows(summary_results)

gene_level_df <- bind_rows(gene_level_results)

# --------------------------------------------------
# Pairwise proportion tests
# --------------------------------------------------

pairwise_prop <- combn(
    seq_len(nrow(summary_df)),
    2,
    simplify = FALSE
)

prop_results <- lapply(pairwise_prop, function(idx) {

    r1 <- summary_df[idx[1],]
    r2 <- summary_df[idx[2],]

    test <- prop.test(
        x = c(
            r1$overlap,
            r2$overlap
        ),
        n = c(
            r1$n_tool_genes,
            r2$n_tool_genes
        )
    )

    data.frame(
        tool_A = r1$tool,
        tool_B = r2$tool,
        p_value = test$p.value
    )
})

prop_results_df <- bind_rows(prop_results)

# --------------------------------------------------
# Logistic regression
# --------------------------------------------------

logit_model <- glm(
    tdp43_target ~ tool,
    data = gene_level_df,
    family = binomial()
)

logit_df <- broom::tidy(
    logit_model,
    exponentiate = TRUE,
    conf.int = TRUE
)

# --------------------------------------------------
# Binding Score
# --------------------------------------------------
# Score = log1p(intron_same_count_density), restricted to genes
# with confirmed intronic TDP-43 crosslinks (density > 0).
# This avoids length confounding (density normalises by segment
# length) and separates binding intensity from binding prevalence.
# Genes with no intronic crosslinks are excluded from the
# intensity comparison; they are captured by the binary
# enrichment analysis (Fisher/logistic) above.

cat(sprintf(
    "\nCalculating binding scores and performing statistical tests...\n"
))

# Build gene -> intronic density lookup.
# Join via SYMBOL only (consistent with iclip_map construction);
# take max density when multiple symbols collapse to one Ensembl ID.
clip_scores <- clip_df %>%
    inner_join(
        iclip_map,
        by = c("gene_name" = "SYMBOL")
    ) %>%
    group_by(ENSEMBL) %>%
    summarise(
        intron_same_count_density = max(intron_same_count_density, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # Retain only genes with genuine intronic crosslink evidence
    filter(
        !is.na(intron_same_count_density),
        intron_same_count_density > 0
    )

# Collect per-tool gene scores (bound genes only)
binding_results <- list()

for (f in tool_files) {

    tool_name <- rename_tool(
        sub("_go_ids\\.txt$", "", basename(f))
    )

    tool_genes <- fread(f)$gene_id %>% unique()

    tool_scores <- data.frame(gene_id = tool_genes) %>%
        inner_join(                          # inner: keep only intronic-bound genes
            clip_scores,
            by = c("gene_id" = "ENSEMBL")
        ) %>%
        mutate(
            binding_score = log1p(intron_same_count_density),
            tool = tool_name
        ) %>%
        dplyr::select(gene_id, tool, intron_same_count_density, binding_score)

    binding_results[[tool_name]] <- tool_scores
}

binding_df <- bind_rows(binding_results)

cat(sprintf(
    "Binding score comparison: %d gene-tool pairs with intron density > 0\n",
    nrow(binding_df)
))

# Global test: do intensity distributions differ across tools?
kruskal_res <- kruskal.test(
    binding_score ~ tool,
    data = binding_df
)

kruskal_df <- data.frame(
    statistic = unname(kruskal_res$statistic),
    df = unname(kruskal_res$parameter),
    p_value = kruskal_res$p.value
)

# Pairwise comparisons (BH-adjusted Wilcoxon tests)
pairwise_res <- pairwise.wilcox.test(
    binding_df$binding_score,
    binding_df$tool,
    p.adjust.method = "BH",
    exact = FALSE
)

pairwise_binding_df <- as.data.frame(as.table(pairwise_res$p.value))
colnames(pairwise_binding_df) <- c("tool_A", "tool_B", "adjusted_p")
pairwise_binding_df <- pairwise_binding_df %>% filter(!is.na(adjusted_p))

# Per-tool summary statistics
binding_summary <- binding_df %>%
    group_by(tool) %>%
    summarise(
        n_intronic_bound = n(),
        median_binding = median(binding_score),
        mean_binding   = mean(binding_score),
        q25 = quantile(binding_score, 0.25),
        q75 = quantile(binding_score, 0.75),
        .groups = "drop"
    )

# --------------------------------------------------
# Visualisation: dotted boxplot of binding scores
# --------------------------------------------------

#tool_order <- c(
#    "W",
#    "W + S",
#    "L",
#    "L + S",
#    "M",
#    "M + S",
#    "Cons",
#    "SpliCeAT"
    #"missed_by_spliceat",    # adjust to whatever rename_tool() produces
    #"all_tools_union"        # adjust likewise
#)
#binding_df$tool <- factor(binding_df$tool, levels = tool_order)

# Annotate n per group for x-axis labels
n_labels <- binding_df %>%
    group_by(tool) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0(tool, "\n(n=", n, ")"))

label_map <- setNames(n_labels$label, n_labels$tool)

p_binding <- ggplot(
    binding_df,
    aes(x = tool, y = binding_score)
) +
    geom_boxplot(
        aes(fill = tool),
        outlier.shape = NA,    # outliers shown via jitter instead
        width = 0.5,
        alpha = 0.6,
        colour = "grey30",
        linewidth = 0.6
    ) +
    geom_jitter(
        width = 0.18,
        size = 0.9,
        alpha = 0.45,
        colour = "grey20"
    ) +
    scale_x_discrete(labels = label_map) +
    scale_fill_brewer(palette = "Set2") +
    labs(
        x        = NULL,
        y        = "TDP-43 intronic binding score\nlog1p(intron crosslink density)",
       # title    = "TDP-43 intronic binding intensity\namong differentially spliced genes"
    ) +
    theme_classic(base_size = 16) +
    theme(
        legend.position    = "none",
        text               = element_text(family = "serif"),
        axis.text.x        = element_text(size = 18, colour = "black", family = "serif"),
        axis.text.y        = element_text(size = 18, colour = "black", family = "serif"),
        axis.title.y       = element_text(size = 20, margin = margin(r = 10), family = "serif"),
       # plot.title         = element_text(size = 16, face = "bold", hjust = 0.5, family = "serif"),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.4)
    )

ggsave(
    file.path(output_dir, "binding_score_boxplot.pdf"),
    plot   = p_binding,
    width  = max(4, length(unique(binding_df$tool)) * 1.4),
    height = 5,
    device = cairo_pdf
)

ggsave(
    file.path(output_dir, "binding_score_boxplot.png"),
    plot  = p_binding,
    width  = max(4, length(unique(binding_df$tool)) * 1.4),
    height = 5,
    dpi   = 300
)

p_violin <- ggplot(
    binding_df,
    aes(x = tool, y = binding_score)
) +
    geom_violin(
        aes(fill = tool),
        alpha      = 0.6,
        colour     = "grey30",
        linewidth  = 0.5,
        trim       = FALSE    # show full distribution tails
    ) +
    geom_boxplot(
        width     = 0.08,     # thin inner boxplot for reference
        outlier.shape = NA,
        colour    = "grey20",
        fill      = "white",
        alpha     = 0.8,
        linewidth = 0.5
    ) +
    scale_x_discrete(labels = label_map) +
    scale_fill_brewer(palette = "Set2") +
    labs(
        x     = NULL,
        y     = "TDP-43 intronic binding score\nlog1p(intron crosslink density)",
       # title = "TDP-43 intronic binding intensity\namong differentially spliced genes"
    ) +
    theme_classic(base_size = 16) +
    theme(
        legend.position    = "none",
        text               = element_text(family = "serif"),
        axis.text.x        = element_text(size = 18, colour = "black", family = "serif"),
        axis.text.y        = element_text(size = 18, colour = "black", family = "serif"),
        axis.title.y       = element_text(size = 20, margin = margin(r = 10), family = "serif"),
       # plot.title         = element_text(size = 16, face = "bold", hjust = 0.5, family = "serif"),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.4)
    )

ggsave(
    file.path(output_dir, "binding_score_violin.pdf"),
    plot   = p_violin,
    width  = max(4, length(unique(binding_df$tool)) * 1.4),
    height = 5,
    device = cairo_pdf
)

ggsave(
    file.path(output_dir, "binding_score_violin.png"),
    plot   = p_violin,
    width  = max(4, length(unique(binding_df$tool)) * 1.4),
    height = 5,
    dpi    = 300
)
# --------------------------------------------------
# Output
# --------------------------------------------------

fwrite(
    summary_df,
    file.path(
        output_dir,
        "tool_enrichment_summary.tsv"
    ),
    sep = "\t"
)

fwrite(
    prop_results_df,
    file.path(
        output_dir,
        "pairwise_proportion_tests.tsv"
    ),
    sep = "\t"
)

fwrite(
    logit_df,
    file.path(
        output_dir,
        "logistic_regression.tsv"
    ),
    sep = "\t"
)

fwrite(
    gene_level_df,
    file.path(
        output_dir,
        "gene_level_table.tsv"
    ),
    sep = "\t"
)

# Binding-score outputs
fwrite(
    binding_df,        # columns: gene_id, tool, intron_same_count_density, binding_score
    file.path(output_dir, "binding_scores.tsv"),
    sep = "\t"
)

fwrite(
    binding_summary,
    file.path(output_dir, "binding_score_summary.tsv"),
    sep = "\t"
)

fwrite(
    pairwise_binding_df,
    file.path(output_dir, "binding_score_pairwise.tsv"),
    sep = "\t"
)

fwrite(
    kruskal_df,
    file.path(output_dir, "binding_score_kruskal.tsv"),
    sep = "\t"
)