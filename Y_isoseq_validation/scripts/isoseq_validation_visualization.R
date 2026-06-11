suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
})
fontfamily <- "Times"
basesize <- 20
fontsize_title <- 26
#'
#' @param validated Data table with validation results
#' @param tool_name Name of the tool
#' @return Data frame with summary statistics
calculate_validation_stats <- function(validated, tool_name) {
  
  overall <- data.frame(
    tool = tool_name,
    event_type = "ALL",
    total = nrow(validated),
    confirmed = sum(validated$isoseq_confirmed),
    pct_confirmed = 100 * sum(validated$isoseq_confirmed) / nrow(validated),
    stringsAsFactors = FALSE
  )
  
  # Per event type
  by_event <- validated %>%
    group_by(event_type) %>%
    summarize(
      total = n(),
      confirmed = sum(isoseq_confirmed),
      pct_confirmed = 100 * sum(isoseq_confirmed) / n(),
      .groups = "drop"
    ) %>%
    mutate(tool = tool_name) %>%
    select(tool, event_type, total, confirmed, pct_confirmed)
  
  return(bind_rows(overall, by_event))
}

#' Print validation summary
#'
#' @param summary_df Summary data frame
print_validation_summary <- function(summary_df) {
  
  cat("\n=== Validation Summary ===\n\n")
  
  # Overall by tool
  overall <- summary_df %>%
    filter(event_type == "ALL") %>%
    arrange(desc(pct_confirmed))
  
  cat("Overall confirmation rates by tool:\n")
  cat(sprintf("%-20s %10s %10s %10s\n", 
              "Tool", "Total", "Confirmed", "% Confirmed"))
  cat(strrep("-", 55), "\n")
  
  for (i in 1:nrow(overall)) {
    cat(sprintf("%-20s %10d %10d %9.1f%%\n",
                overall$tool[i],
                overall$total[i],
                overall$confirmed[i],
                overall$pct_confirmed[i]))
  }
  
  cat("\n")
  
  # By event type (aggregate across tools)
  by_event <- summary_df %>%
    filter(event_type != "ALL") %>%
    group_by(event_type) %>%
    summarize(
      total = sum(total),
      confirmed = sum(confirmed),
      pct_confirmed = 100 * sum(confirmed) / sum(total)
    ) %>%
    arrange(desc(pct_confirmed))
  
  cat("Confirmation rates by event type (all tools):\n")
  cat(sprintf("%-35s %10s %10s %10s\n", 
              "Event Type", "Total", "Confirmed", "% Confirmed"))
  cat(strrep("-", 70), "\n")
  
  for (i in 1:nrow(by_event)) {
    cat(sprintf("%-35s %10d %10d %9.1f%%\n",
                by_event$event_type[i],
                by_event$total[i],
                by_event$confirmed[i],
                by_event$pct_confirmed[i]))
  }
  
  cat("\n")
}

#' Generate validation report plots
#'
#' @param summary_df Summary data frame
#' @param output_dir Output directory
#' Generate validation report plots
#'
#' @param summary_df Summary data frame
#' @param output_dir Output directory
generate_validation_plots <- function(summary_df, output_dir = ".") {
  # 1. Rename tools for visualization
  summary_df <- summary_df %>%
    mutate(tool = case_when(
      tool == "Majiq_validated" ~ "MAJIQ",
      tool == "Whippet_validated" ~ "Whippet",
      tool == "LeafCutter_validated" ~ "LeafCutter",
      tool == "Majiq_+_StringTie_validated" ~ "MAJIQ + StringTie",
      tool == "Whippet_+_StringTie_validated" ~ "Whippet + StringTie",
      tool == "LeafCutter_+_StringTie_validated" ~ "LeafCutter + StringTie",
      tool == "DS_consensus_validated" ~ "DS consensus",
      tool == "SpliCeAT_validated" ~ "SpliCeAT",
      TRUE ~ tool
    ))

  # 2. Define tool order for plotting
  tool_order <- c("SpliCeAT", "DS consensus", "MAJIQ + StringTie",
                  "MAJIQ", "LeafCutter + StringTie", "LeafCutter",
                  "Whippet + StringTie", "Whippet")

  # 3. FIX: Convert 'tool' to a factor with explicitly defined levels
  summary_df <- summary_df %>%
    mutate(tool = factor(tool, levels = tool_order))

  library(ggplot2)
  
  # Plot 1: Overall confirmation by tool
  overall <- summary_df %>% filter(event_type == "ALL")
  
  p1 <- ggplot(overall, aes(x = reorder(tool, pct_confirmed), 
                           y = pct_confirmed, fill = tool)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", pct_confirmed)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    ylim(0, 100) +
    theme_minimal() +
    labs(title = "IsoSeq Confirmation Rate by Tool",
         x = "Tool", y = "% Confirmed") +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "validation_by_tool.pdf"), 
         p1, width = 8, height = 6)
  
  # Plot 2: Confirmation by event type
  by_event <- summary_df %>% 
    filter(event_type != "ALL") %>%
    group_by(event_type) %>%
    summarize(
      total = sum(total),
      confirmed = sum(confirmed),
      pct_confirmed = 100 * sum(confirmed) / sum(total)
    )
  
  p2 <- ggplot(by_event, aes(x = reorder(event_type, pct_confirmed),
                            y = pct_confirmed)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", pct_confirmed)),
              hjust = -0.2, size = 3) +
    coord_flip() +
    ylim(0, 100) +
    theme_minimal() +
    labs(title = "IsoSeq Confirmation Rate by Event Type",
         x = "Event Type", y = "% Confirmed")
  
  ggsave(file.path(output_dir, "validation_by_event.pdf"),
         p2, width = 10, height = 6)
  
  # Plot 3: Heatmap of tool vs event type
  heatmap_data <- summary_df %>%
    filter(event_type != "ALL") %>%
    select(tool, event_type, pct_confirmed)
  
  p3 <- ggplot(heatmap_data, aes(x = event_type, y = tool, 
                                fill = pct_confirmed)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.0f%%", pct_confirmed)), size = 3) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                        midpoint = 50, limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "IsoSeq Confirmation: Tool vs Event Type",
         x = "Event Type", y = "Tool", fill = "% Confirmed")
  
  ggsave(file.path(output_dir, "validation_heatmap.pdf"),
         p3, width = 12, height = 6)
  
# Plot 4: Concentric bubble chart of total vs confirmed events
bubble_data <- summary_df %>%
  filter(event_type != "ALL")

  p4 <- ggplot(bubble_data, aes(x = event_type, y = tool)) +
    # Outer circle (Total)
    geom_point(aes(size = total), shape = 21, color = "grey40", fill = NA) +
    
    # Inner circle (Confirmed) 
    geom_point(aes(size = confirmed, color = pct_confirmed), shape = 16) +
    
    # Color Gradient (will now apply correctly to the inner circle)
    scale_color_gradient2(
      low = "purple", mid = "grey90", high = "darkgreen",
      midpoint = 50, limits = c(0, 100),
      name = "% Confirmed" # Adds a nice title to your color legend
    ) +
    scale_size_area(max_size = 12) +
    theme_minimal() +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y       = element_text(color = "black"),
      panel.grid.major  = element_line(color = "grey95"),
      text              = element_text(family = fontfamily, size = basesize)
    ) +
    labs(
      x = "Event Type",
      y = "Tool",
      size = "Number of Events"
    )

  ggsave(file.path(output_dir, "validation_concentric_bubbles.pdf"),
        p4, width = 12, height = 6)
    
    cat(sprintf("Plots saved to %s\n", output_dir))
  }

  #' Collapse cEI clusters into single events
  #'
  #' For each cluster where cEI == TRUE, create a single representative event
  #' with event_type = "cEI". The collapsed event is confirmed only if ALL
  #' original events in the cluster were confirmed.
  #'
  #' @param validated Data table with validation results (must have cluster_id and cEI columns)
  #' @return Data table with cEI clusters collapsed
  collapse_cEI_clusters <- function(validated) {
    
    if (!"cEI" %in% colnames(validated) || !"cluster_id" %in% colnames(validated)) {
      return(validated)
    }
    
    # Separate cEI and non-cEI events
    non_cei_events <- validated[is.na(validated$cEI) | validated$cEI == FALSE]
    cei_clusters <- validated[validated$cEI == TRUE]
    
    if (nrow(cei_clusters) == 0) {
      return(validated)
    }
    
    # Group cEI events by cluster_id
    cei_clusters_by_id <- cei_clusters %>%
      group_by(cluster_id) %>%
      group_split()
    
    # Collapse each cEI cluster to a single event
    collapsed_events <- list()
    
    for (cluster_group in cei_clusters_by_id) {
      cluster_id <- cluster_group$cluster_id[1]
      
      # Check if all events in this cluster are confirmed
      all_confirmed <- all(cluster_group$isoseq_confirmed == TRUE, na.rm = TRUE)
      
      # Use first event as template and update fields
      collapsed_event <- cluster_group[1, ]
      collapsed_event$event_type <- "cEI"
      collapsed_event$isoseq_confirmed <- all_confirmed
      
      # Combine transcript IDs from all events in the cluster
      all_transcripts <- unique(unlist(strsplit(
        paste(cluster_group$isoseq_transcript_ids, collapse = ","), ","
      )))
      all_transcripts <- all_transcripts[!is.na(all_transcripts) & nchar(all_transcripts) > 0]
      collapsed_event$isoseq_transcript_ids <- paste(all_transcripts, collapse = ",")
      
      collapsed_events[[length(collapsed_events) + 1]] <- collapsed_event
    }
    
    # Combine collapsed cEI events with non-cEI events
    if (length(collapsed_events) > 0) {
      collapsed_df <- rbindlist(collapsed_events)
      result <- rbind(non_cei_events, collapsed_df)
    } else {
      result <- non_cei_events
    }
    
    cat(sprintf("      Collapsed %d cEI clusters to %d representative events\n",
                nrow(cei_clusters), length(collapsed_events)))
    
    return(result)
  }

# ============================================================================
# Main Snakemake Execution Block
# ============================================================================

tryCatch({
  
  # Get input/output from Snakemake
  input_files <- snakemake@input
  output_file <- snakemake@output[[1]]
  output_dir <- snakemake@params[["output_dir"]]
  
  # Optional: get log file if available
  log_file <- NULL
  if (!is.null(snakemake@log) && length(snakemake@log) > 0) {
    log_file <- snakemake@log[[1]]
    log_con <- file(log_file, open = "wt")
    sink(log_con, split = TRUE)
    sink(log_con, type = "message")
  }
  
  cat("=== IsoSeq Validation Visualization ===\n\n")
  cat(sprintf("Loading %d validation output files...\n", length(input_files)))
  cat(sprintf("Output directory: %s\n\n", output_dir))
  
  # Create output directory if needed
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load all validation files
  all_validated <- list()
  all_stats <- list()
  
  for (i in seq_along(input_files)) {
    file <- input_files[[i]]
    tool_name <- tools::file_path_sans_ext(basename(file))
    
    cat(sprintf("[%2d/%d] Loading %s...\n", i, length(input_files), tool_name))
    
    tryCatch({
      validated <- fread(file)
      
      # Collapse cEI clusters if present
      if ("cEI" %in% colnames(validated) && "cluster_id" %in% colnames(validated)) {
        cat(sprintf("      Collapsing cEI clusters...\n"))
        validated <- collapse_cEI_clusters(validated)
      }
      
      all_validated[[tool_name]] <- validated
      
      # Calculate stats for this tool
      stats <- calculate_validation_stats(validated, tool_name)
      all_stats[[tool_name]] <- stats
      
      # Print tool summary
      n_confirmed <- sum(validated$isoseq_confirmed == TRUE, na.rm = TRUE)
      pct_confirmed <- 100 * n_confirmed / nrow(validated)
      cat(sprintf("      ✓ %3d / %-3d confirmed (%.1f%%) \n", 
                  n_confirmed, nrow(validated), pct_confirmed))
      
    }, error = function(e) {
      cat(sprintf("      ✗ ERROR: %s\n", conditionMessage(e)))
    })
  }
  
  if (length(all_stats) == 0) {
    stop("No validation files were successfully loaded")
  }
  
  cat("\n")
  
  # Combine statistics
  combined_stats <- bind_rows(all_stats)
  
  # Print summary
  print_validation_summary(combined_stats)
  
  # Generate plots
  cat("Generating visualization plots...\n")
  generate_validation_plots(combined_stats, output_dir)
  
  # Write comprehensive summary statistics
  summary_file <- file.path(output_dir, "validation_summary_statistics.tsv")
  cat(sprintf("Writing summary statistics to: %s\n", summary_file))
  fwrite(combined_stats, summary_file, sep = "\t", quote = FALSE)
  
  # Write detailed tool comparison
  cat("\nGenerating detailed comparison files...\n\n")
  
  # 1. Per-tool totals
  tool_totals <- combined_stats %>%
    filter(event_type == "ALL") %>%
    select(tool, total, confirmed, pct_confirmed) %>%
    arrange(desc(pct_confirmed))
  
  tool_summary_file <- file.path(output_dir, "tool_summary.tsv")
  fwrite(tool_totals, tool_summary_file, sep = "\t", quote = FALSE)
  cat(sprintf("Saved tool summary: %s\n", tool_summary_file))
  
  # 2. Event type comparison
  event_totals <- combined_stats %>%
    filter(event_type != "ALL") %>%
    group_by(event_type) %>%
    summarize(
      total_events = sum(total),
      total_confirmed = sum(confirmed),
      avg_pct_confirmed = mean(pct_confirmed)
    ) %>%
    arrange(desc(total_confirmed))
  
  event_summary_file <- file.path(output_dir, "event_type_summary.tsv")
  fwrite(event_totals, event_summary_file, sep = "\t", quote = FALSE)
  cat(sprintf("Saved event type summary: %s\n", event_summary_file))
  
  # 3. Tool-by-event-type matrix
  event_matrix <- combined_stats %>%
    filter(event_type != "ALL") %>%
    select(tool, event_type, pct_confirmed) %>%
    pivot_wider(names_from = event_type, values_from = pct_confirmed)
  
  matrix_file <- file.path(output_dir, "tool_by_event_matrix.tsv")
  fwrite(event_matrix, matrix_file, sep = "\t", quote = FALSE)
  cat(sprintf("Saved tool-by-event matrix: %s\n", matrix_file))
  
  # Write main summary output file
  cat(sprintf("\nWriting main summary output to: %s\n", output_file))
  summary_text <- paste0(
    "=== IsoSeq Validation Summary Report ===\n\n",
    
    "OVERALL STATISTICS:\n",
    sprintf("Number of tools validated: %s\n", length(all_validated)),
    sprintf("Total events across all tools: %d\n", sum(tool_totals$total)),
    sprintf("Total confirmed events: %d\n", sum(tool_totals$confirmed)),
    sprintf("Overall confirmation rate: %.1f%%\n\n",
            100 * sum(tool_totals$confirmed) / sum(tool_totals$total)),
    
    "TOOL PERFORMANCE:\n",
    paste(apply(tool_totals, 1, function(row) {
      sprintf("  %s: %d/%d (%.1f%%)",
        row["tool"],
        as.integer(row["confirmed"]),
        as.integer(row["total"]),
        as.numeric(row["pct_confirmed"]))
    }), collapse = "\n"),
    "\n\n",
    
    "EVENT TYPE PERFORMANCE:\n",
    paste(apply(event_totals, 1, function(row) {
      sprintf("  %s: %d confirmed out of %d events (%.1f%% avg)",
        row["event_type"],
        as.integer(row["total_confirmed"]),
        as.integer(row["total_events"]),
        as.numeric(row["avg_pct_confirmed"]))
    }), collapse = "\n"),
    "\n\n",
    
    "OUTPUT FILES:\n",
    sprintf("  Summary statistics: %s\n", summary_file),
    sprintf("  Tool summary: %s\n", tool_summary_file),
    sprintf("  Event type analysis: %s\n", event_summary_file),
    sprintf("  Tool-by-event matrix: %s\n", matrix_file),
    sprintf("  Visualizations: %s/*.pdf\n", output_dir),
    sprintf("  Log file: %s\n\n", log_file),
    
    "RECOMMENDATIONS:\n",
    if_else(mean(tool_totals$pct_confirmed) < 10,
            "  ⚠ Low overall confirmation rate - review alignment tolerance\n",
            "  ✓ Reasonable confirmation rates\n"),
    if_else(max(tool_totals$pct_confirmed) > 2 * min(tool_totals$pct_confirmed),
            "  ⚠ Large variation between tools - check for tool-specific biases\n",
            "  ✓ Consistent performance across tools\n")
  )
  
  writeLines(summary_text, output_file)
  
  cat("\n✓ Visualization complete!\n")
  cat(sprintf("Summary report written to: %s\n", output_file))
  
  # Close log file
  if (!is.null(log_file)) {
    sink()
    sink(type = "message")
    close(log_con)
  }
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat(traceback(), "\n")
  if (!is.null(log_file)) {
    sink()
    sink(type = "message")
    close(log_con)
  }
  quit(status = 1)
})

