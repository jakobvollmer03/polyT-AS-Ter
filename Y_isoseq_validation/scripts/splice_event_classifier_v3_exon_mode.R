#!/usr/bin/env Rscript
## This script classifies exon contigs (partial/full exons) produced by a splice
## detection tool.  Unlike the ground-truth version this script accepts a SINGLE
## tool-output file that contains a 'dpsi' column:
##
##   dpsi > 0  →  contig is sample-specific          (equivalent to "up")
##   dpsi < 0  →  contig is reference-specific       (equivalent to "down")
##
## Classification rules:
##
## "up" contigs (dpsi > 0):
##   - spans intronic region, connects to exon on 5' end → alternative_donor_late
##   - spans intronic region, connects to exon on 3' end → alternative_acceptor_early
##   - spans intronic region, connects to exons on both ends → intron_retention
##   - spans intronic region, no connections, within known gene span → cryptic_exon
##   - spans intronic region, no connections, outside any gene span → novel_contig (discarded)
##   - spans exonic region → ERROR (discarded)
##
## "down" contigs (dpsi < 0):
##   - spans exonic region, connects to 5' end of exon → alternative_acceptor_late
##   - spans exonic region, connects to 3' end of exon → alternative_donor_early
##   - spans exonic region, connects to neither boundary → exonic_cryptic
##   - spans intronic region → ERROR (discarded)
##
## Both ERROR_* and novel_contig events are discarded before output is written.
##
## Input table format (tool output): chr, start, end, strand, dpsi  [+ optional extra columns]
## Snakemake params mirror splice_event_classifier_v3.R:
##   params: junction_file (or novel_features), reference_gtf (or gtf_file), tolerance
##   output: [1]  classified table
##   log:    [1]  log file

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Snakemake I/O  — mirrors splice_event_classifier_v3.R exactly
# ---------------------------------------------------------------------------
if (!is.null(snakemake@params[["junction_file"]])) {
  junction_file <- snakemake@params[["junction_file"]]
} else if (!is.null(snakemake@params[["novel_features"]])) {
  junction_file <- snakemake@params[["novel_features"]]
} else {
  junction_file <- snakemake@params[[1]]   # positional fallback
}

gtf_file <- if (!is.null(snakemake@params[["gtf_file"]])) {
  snakemake@params[["gtf_file"]]
} else if (!is.null(snakemake@params[["reference_gtf"]])) {
  snakemake@params[["reference_gtf"]]
} else {
  snakemake@params[[2]]                    # positional fallback
}

output_file <- snakemake@output[[1]]
log_file    <- snakemake@log[[1]]
tolerance   <- if (!is.null(snakemake@params[["tolerance"]])) {
  as.numeric(snakemake@params[["tolerance"]])
} else {
  5
}

# ---------------------------------------------------------------------------
# Redirect stdout + stderr to log file
# ---------------------------------------------------------------------------
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}

cat("=== Exon Contig Classifier (tool-output mode) ===\n\n")
cat(sprintf("Input file  : %s\n", junction_file))
cat(sprintf("GTF file    : %s\n", gtf_file))
cat(sprintf("Output file : %s\n", output_file))
cat(sprintf("Tolerance   : %d bp\n\n", tolerance))

# ===========================================================================
# FUNCTION DEFINITIONS
# ===========================================================================

#' Load tool-output exon contig file and split on dpsi sign.
#'
#' The file must contain at minimum: chr, start, end, strand, dpsi.
#' Rows with dpsi > 0 are treated as "up" (sample-specific).
#' Rows with dpsi < 0 are treated as "down" (reference-specific).
#' Rows with dpsi == 0 or NA are skipped with a warning.
#'
#' @param contig_file Path to tool-output file (TSV, header required).
#' @return Named list with elements $up and $down, each a GRanges.
load_tool_output <- function(contig_file) {
  cat(sprintf("Loading tool output from: %s\n", contig_file))

  if (!file.exists(contig_file)) {
    stop(sprintf("Input file not found: %s", contig_file))
  }

  df <- as.data.table(read.table(contig_file,
                                  header           = TRUE,
                                  stringsAsFactors = FALSE,
                                  comment.char     = "#",
                                  fill             = TRUE,
                                  row.names        = NULL))

  required_cols <- c("chr", "start", "end", "strand", "dpsi")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Tool output is missing required column(s): %s",
                 paste(missing, collapse = ", ")))
  }

  df$dpsi  <- as.numeric(df$dpsi)
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)

  n_zero <- sum(df$dpsi == 0 | is.na(df$dpsi), na.rm = TRUE)
  if (n_zero > 0) {
    warning(sprintf("%d row(s) with dpsi == 0 or NA will be skipped.", n_zero))
    df <- df[!is.na(df$dpsi) & df$dpsi != 0, ]
  }

  cat(sprintf("  Total rows loaded: %d\n", nrow(df)))

  df_up   <- df[df$dpsi > 0, ]
  df_down <- df[df$dpsi < 0, ]

  cat(sprintf("  dpsi > 0 (up)   : %d rows\n", nrow(df_up)))
  cat(sprintf("  dpsi < 0 (down) : %d rows\n", nrow(df_down)))

  make_gr <- function(d) {
    if (nrow(d) == 0) return(GRanges())
    gr <- GRanges(
      seqnames = d$chr,
      ranges   = IRanges(start = d$start, end = d$end),
      strand   = d$strand
    )
    # Carry all extra metadata columns (e.g. dpsi, gene_id from tool) through
    extra_cols <- setdiff(colnames(d), c("chr", "start", "end", "strand"))
    if (length(extra_cols) > 0) {
      mcols(gr) <- d[, extra_cols, with = FALSE]
    }
    gr
  }

  list(up   = make_gr(df_up),
       down = make_gr(df_down))
}

#' Load reference annotation
load_reference <- function(gtf_file) {
  cat("Loading reference annotation...\n")

  gtf   <- rtracklayer::import(gtf_file)
  exons <- gtf[gtf$type == "exon"]

  cat(sprintf("Loaded %d exons from %d transcripts\n",
              length(exons),
              length(unique(exons$transcript_id))))

  return(exons)
}

#' Determine whether a contig region is exonic in the reference.
determine_region_type <- function(contig, exons_gr) {

  chr           <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))

  chr_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                        as.character(strand(exons_gr)) == contig_strand]

  if (length(chr_exons) == 0) {
    return(list(is_exonic = FALSE, overlapping_exons = GRanges()))
  }

  overlapping <- chr_exons[overlapsAny(chr_exons, contig)]

  if (length(overlapping) == 0) {
    return(list(is_exonic = FALSE, overlapping_exons = GRanges()))
  }

  exon_union     <- reduce(overlapping)
  contig_covered <- all(overlapsAny(contig, exon_union, type = "within"))

  return(list(is_exonic = contig_covered, overlapping_exons = overlapping))
}

#' Check whether a contig falls within the genomic span of any known gene.
#'
#' A gene's span is defined as [min(exon_start), max(exon_end)] across all of
#' its exons on the same chr/strand.  This is used to distinguish intronic
#' cryptic exons (within a gene) from completely intergenic novel contigs.
#'
#' @param contig    Single-element GRanges.
#' @param exons_gr  Reference exons GRanges (all chromosomes).
#' @return Named list: $within_gene (logical), $gene_id (character or NA).
is_within_known_gene <- function(contig, exons_gr) {

  chr           <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  contig_start  <- start(contig)
  contig_end    <- end(contig)

  chr_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                        as.character(strand(exons_gr)) == contig_strand]

  if (length(chr_exons) == 0) {
    return(list(within_gene = FALSE, gene_id = NA_character_))
  }

  # Build a per-gene span table
  gene_ids     <- chr_exons$gene_id
  unique_genes <- unique(gene_ids)

  for (gid in unique_genes) {
    g_exons    <- chr_exons[chr_exons$gene_id == gid]
    gene_start <- min(start(g_exons))
    gene_end   <- max(end(g_exons))
    if (contig_start >= gene_start && contig_end <= gene_end) {
      return(list(within_gene = TRUE, gene_id = gid))
    }
  }

  return(list(within_gene = FALSE, gene_id = NA_character_))
}

#' Find the exon the contig connects to at its biological 5' end.
find_upstream_connection <- function(contig, exons_gr, tolerance = 5) {

  chr           <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  contig_start  <- start(contig)
  contig_end    <- end(contig)

  gene_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                         as.character(strand(exons_gr)) == contig_strand]

  if (contig_strand == "+") {
    nearby <- gene_exons[end(gene_exons) < contig_start &
                         end(gene_exons) >= (contig_start - tolerance)]
  } else {
    nearby <- gene_exons[start(gene_exons) > contig_end &
                         start(gene_exons) <= (contig_end + tolerance)]
  }

  if (length(nearby) == 0) return(data.table())

  closest_idx <- if (contig_strand == "+") which.max(end(nearby)) else which.min(start(nearby))
  closest     <- nearby[closest_idx]

  data.table(
    gene_id       = closest$gene_id,
    transcript_id = closest$transcript_id,
    exon_boundary = if (contig_strand == "+") end(closest) else start(closest),
    distance      = if (contig_strand == "+") {
      contig_start - end(closest)
    } else {
      start(closest) - contig_end
    }
  )
}

#' Find the exon the contig connects to at its biological 3' end.
find_downstream_connection <- function(contig, exons_gr, tolerance = 5) {

  chr           <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  contig_start  <- start(contig)
  contig_end    <- end(contig)

  gene_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                         as.character(strand(exons_gr)) == contig_strand]

  if (contig_strand == "+") {
    nearby <- gene_exons[start(gene_exons) > contig_end &
                         start(gene_exons) <= (contig_end + tolerance)]
  } else {
    nearby <- gene_exons[end(gene_exons) < contig_start &
                         end(gene_exons) >= (contig_start - tolerance)]
  }

  if (length(nearby) == 0) return(data.table())

  closest_idx <- if (contig_strand == "+") which.min(start(nearby)) else which.max(end(nearby))
  closest     <- nearby[closest_idx]

  data.table(
    gene_id       = closest$gene_id,
    transcript_id  = closest$transcript_id,
    exon_boundary  = if (contig_strand == "+") start(closest) else end(closest),
    distance       = if (contig_strand == "+") {
      start(closest) - contig_end
    } else {
      contig_start - end(closest)
    }
  )
}

#' Classify a single exon contig.
#'
#' @param contig_idx  Integer index into contig_gr.
#' @param contig_gr   GRanges with all contigs of one type.
#' @param exons_gr    GRanges with reference exons.
#' @param contig_type "up" (dpsi > 0) or "down" (dpsi < 0).
#' @param tolerance   Coordinate-matching tolerance in bp.
#' @return Single-row data.frame with classification.
classify_exon_contig <- function(contig_idx, contig_gr, exons_gr,
                                 contig_type = "up", tolerance = 5) {

  contig <- contig_gr[contig_idx]

  result <- data.frame(
    contig_id     = sprintf("%s:%d-%d:%s_%d",
                            seqnames(contig), start(contig), end(contig),
                            strand(contig), contig_idx),
    chr           = as.character(seqnames(contig)),
    start         = start(contig),
    end           = end(contig),
    strand        = as.character(strand(contig)),
    contig_type   = contig_type,
    event_type    = "unclassified",
    gene_id       = NA_character_,
    transcript_id = NA_character_,
    confidence    = NA_character_,
    stringsAsFactors = FALSE
  )

  # Carry dpsi through if present in the GRanges metadata
  if ("dpsi" %in% colnames(mcols(contig_gr))) {
    result$dpsi <- mcols(contig_gr)$dpsi[contig_idx]
  }

  region_info       <- determine_region_type(contig, exons_gr)
  is_exonic         <- region_info$is_exonic
  overlapping_exons <- region_info$overlapping_exons

  # ---- "up" contigs (dpsi > 0, sample-specific) ---------------------------
  if (contig_type == "up") {

    # Up contigs must be intronic; flag if they overlap an annotated exon
    if (is_exonic && length(overlapping_exons) > 0) {
      result$event_type <- "ERROR_spans_exonic"
      result$confidence <- "0"
      return(result)
    }

    upstream_conn   <- find_upstream_connection(contig, exons_gr, tolerance)
    downstream_conn <- find_downstream_connection(contig, exons_gr, tolerance)

    has_upstream   <- nrow(upstream_conn) > 0
    has_downstream <- nrow(downstream_conn) > 0

    if (has_upstream && !has_downstream) {
      result$event_type    <- "cADds"
      result$gene_id       <- upstream_conn$gene_id[1]
      result$transcript_id <- upstream_conn$transcript_id[1]
      result$confidence    <- "1"
    } else if (!has_upstream && has_downstream) {
      result$event_type    <- "cAAus"
      result$gene_id       <- downstream_conn$gene_id[1]
      result$transcript_id <- downstream_conn$transcript_id[1]
      result$confidence    <- "1"
    } else if (has_upstream && has_downstream) {
      result$event_type <- "cIR"
      result$gene_id    <- upstream_conn$gene_id[1]
      result$transcript_id <- upstream_conn$transcript_id[1]
      result$confidence <- "1"
    } else {
      # No connection to any annotated exon boundary.
      # Only assign cryptic_exon if the contig falls within the genomic span
      # (min exon start → max exon end) of a known gene on the same chr/strand.
      # Contigs outside any gene span lack the necessary context and are
      # labelled novel_contig so they can be discarded cleanly.
      gene_check <- is_within_known_gene(contig, exons_gr)
      if (gene_check$within_gene) {
        result$event_type <- "cEI"
        result$gene_id    <- gene_check$gene_id
        result$confidence <- "2"
      } else {
        result$event_type <- "novel_contig"
        result$confidence <- "3"
      }
    }

  # ---- "down" contigs (dpsi < 0, reference-specific) ----------------------
  } else if (contig_type == "down") {

    # Down contigs must be exonic; flag if they fall entirely in an intron
    if (!is_exonic && length(overlapping_exons) == 0) {
      result$event_type <- "ERROR_spans_intronic"
      result$confidence <- "0"
      return(result)
    }

    if (length(overlapping_exons) == 0) {
      result$event_type <- "unclassified"
      result$confidence <- "0"
      return(result)
    }

    ref_exon      <- overlapping_exons[1]
    contig_strand <- as.character(strand(contig))
    contig_start  <- start(contig)
    contig_end    <- end(contig)
    ref_start     <- start(ref_exon)
    ref_end       <- end(ref_exon)

    # Strand-aware 5'/3' boundary matching:
    # + strand: biological 5' = ref_start (acceptor), biological 3' = ref_end (donor)
    # - strand: biological 5' = ref_end   (acceptor), biological 3' = ref_start (donor)
    if (contig_strand == "+") {
      five_prime_match  <- abs(contig_start - ref_start) <= tolerance
      three_prime_match <- abs(contig_end   - ref_end)   <= tolerance
    } else {
      five_prime_match  <- abs(contig_end   - ref_end)   <= tolerance
      three_prime_match <- abs(contig_start - ref_start) <= tolerance
    }

    if (five_prime_match && !three_prime_match) {
      result$event_type    <- "cAAds"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "1"
    } else if (!five_prime_match && three_prime_match) {
      result$event_type    <- "cADus"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "1"
    } else if (!five_prime_match && !three_prime_match) {
      result$event_type    <- "cIJ"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "2"
    } else {
      # Both boundaries match → exon-skipping placeholder
      result$event_type    <- "cES1"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "3"
    }
  }

  return(result)
}

# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

tryCatch({

  cat("\n=== Loading Data ===\n")

  exons_gr <- load_reference(gtf_file)
  cat("\n")

  contigs          <- load_tool_output(junction_file)
  all_contigs_up   <- contigs$up
  all_contigs_down <- contigs$down

  # ---- Classify contigs ----------------------------------------------------

  cat("\n=== Classifying Exon Contigs ===\n")

  results    <- list()
  result_idx <- 1

  if (length(all_contigs_up) > 0) {
    cat(sprintf("\nClassifying %d up contigs (dpsi > 0)...\n", length(all_contigs_up)))
    pb <- txtProgressBar(min = 0, max = length(all_contigs_up), style = 3)
    for (i in seq_along(all_contigs_up)) {
      results[[result_idx]] <- classify_exon_contig(i, all_contigs_up, exons_gr,
                                                    "up", tolerance)
      result_idx <- result_idx + 1
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  if (length(all_contigs_down) > 0) {
    cat(sprintf("\nClassifying %d down contigs (dpsi < 0)...\n", length(all_contigs_down)))
    pb <- txtProgressBar(min = 0, max = length(all_contigs_down), style = 3)
    for (i in seq_along(all_contigs_down)) {
      results[[result_idx]] <- classify_exon_contig(i, all_contigs_down, exons_gr,
                                                    "down", tolerance)
      result_idx <- result_idx + 1
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  # ---- Combine results -----------------------------------------------------

  cat("\n\n=== Combining Results ===\n")

  if (length(results) == 0) {
    cat("No contigs to classify!\n")
    results_df <- data.frame()
  } else {
    results_df <- bind_rows(results)
  }

  # ---- Discard ERROR events ------------------------------------------------
  # ERROR_spans_exonic  : up contig that overlaps an annotated exon (invalid)
  # ERROR_spans_intronic: down contig that falls entirely in an intron (invalid)
  error_types <- c("ERROR_spans_exonic", "ERROR_spans_intronic")
  n_error <- sum(results_df$event_type %in% error_types, na.rm = TRUE)
  if (n_error > 0) {
    cat(sprintf("\nDiscarding %d ERROR event(s) (correspond to canonical exonic/intronic sequences).\n",
                n_error))
    results_df <- results_df[!(results_df$event_type %in% error_types), ]
  }

  # ---- Discard novel contigs -----------------------------------------------
  # novel_contig: up contig with no exon connections AND outside any gene span.
  # These lack the genomic context required for meaningful classification.
  n_novel <- sum(results_df$event_type == "novel_contig", na.rm = TRUE)
  if (n_novel > 0) {
    cat(sprintf("\nDiscarding %d novel_contig(s) with no nearby annotation.\n", n_novel))
    results_df <- results_df[!(results_df$event_type == "novel_contig"), ]
  }

  # ---- Summary statistics --------------------------------------------------

  cat("\n=== Classification Summary ===\n")

  if (nrow(results_df) > 0) {
    print(table(results_df$event_type))

    cat("\nBreakdown by contig type (up/down):\n")
    print(table(results_df$contig_type, results_df$event_type))
  } else {
    cat("No results to summarize\n")
  }

  # ---- Write output --------------------------------------------------------

  cat(sprintf("\nWriting results to %s\n", output_file))
  fwrite(results_df, output_file, sep = "\t", quote = FALSE)

  cat("\n✓ Exon contig classification complete!\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Close log file
if (!is.null(log_file)) {
  sink()
  sink(type = "message")
  close(log_con)
}