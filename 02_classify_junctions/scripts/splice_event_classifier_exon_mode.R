#!/usr/bin/env Rscript
## This script classifies exon contigs (partial/full exons unique to sample or reference).
## The classification differentiates between exon_contigs_up (sample-specific) and 
## exon_contigs_down (reference-specific).
##
## exon_contigs_up classification:
##   - spans intronic region, connects to exon on 5' end → alternative_donor_late
##   - spans intronic region, connects to exon on 3' end → alternative_acceptor_early
##   - spans intronic region, connects to no exon → cryptic_exon
##   - spans exonic region → ERROR
##
## exon_contigs_down classification:
##   - spans exonic region, connects to 5' end of exon → alternative_acceptor_late
##   - spans exonic region, connects to 3' end of exon → alternative_donor_early
##   - spans exonic region, connects to none of exon boundaries → exonic_cryptic
##   - spans intronic region → ERROR

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(data.table)
})

# Get input/output from Snakemake
exon_contigs_up_file <- snakemake@input[["exon_contigs_up"]]
exon_contigs_down_file <- snakemake@input[["exon_contigs_down"]]
gtf_file <- snakemake@params[["reference_gtf"]]
output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]
tolerance <- if (!is.null(snakemake@params[["tolerance"]])) {
  as.numeric(snakemake@params[["tolerance"]])
} else {
  2
}

# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}

cat("=== Exon Contig Classifier ===\n\n")
cat(sprintf("Tolerance for coordinate matching: %d bp\n\n", tolerance))

# ===== FUNCTION DEFINITIONS =====

#' Load exon contigs from file
#'
#' @param contig_file Path to exon contig file
#' @param contig_type "up" or "down"
#' @return GRanges object with contigs
load_exon_contigs <- function(contig_file, contig_type = "up") {
  cat(sprintf("Loading exon contigs (%s) from: %s\n", contig_type, contig_file))
  
  if (!file.exists(contig_file)) {
    cat(sprintf("Warning: Exon contig file not found: %s\n", contig_file))
    return(GRanges())
  }
  
  # Read file
  contigs_df <- read.table(contig_file, 
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           comment.char = "#")
  
  if (nrow(contigs_df) == 0) {
    cat(sprintf("No contigs in %s file\n", contig_type))
    return(GRanges())
  }
  
  # Convert to GRanges
  contigs_gr <- GRanges(
    seqnames = contigs_df$chr,
    ranges = IRanges(start = contigs_df$start, end = contigs_df$end),
    strand = contigs_df$strand
  )
  
  cat(sprintf("Loaded %d contigs (%s)\n", length(contigs_gr), contig_type))
  
  return(contigs_gr)
}

#' Load reference annotation
#'
#' @param gtf_file Path to GTF file
#' @return GRanges object with exons
load_reference <- function(gtf_file) {
  cat("Loading reference annotation...\n")
  
  gtf <- rtracklayer::import(gtf_file)
  exons <- gtf[gtf$type == "exon"]
  
  cat(sprintf("Loaded %d exons from %d transcripts\n", 
              length(exons),
              length(unique(exons$transcript_id))))
  
  return(exons)
}

#' Check if a contig region is exonic or intronic in reference
#'
#' @param contig GRanges object with single contig
#' @param exons_gr GRanges object with reference exons
#' @return List with elements: is_exonic (logical), overlapping_exons (GRanges)
determine_region_type <- function(contig, exons_gr) {
  
  chr <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  
  # Pre-filter to chr/strand before any overlap operation to avoid seqlevel mismatch
  chr_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                        as.character(strand(exons_gr)) == contig_strand]
  
  if (length(chr_exons) == 0) {
    return(list(is_exonic = FALSE, overlapping_exons = GRanges()))
  }
  
  overlapping <- chr_exons[overlapsAny(chr_exons, contig)]
  
  if (length(overlapping) == 0) {
    return(list(is_exonic = FALSE, overlapping_exons = GRanges()))
  }
  
  exon_union <- reduce(overlapping)
  contig_covered <- all(overlapsAny(contig, exon_union, type = "within"))
  
  return(list(is_exonic = contig_covered, overlapping_exons = overlapping))
}

#' Find what exon the contig connects to at its 5' end (upstream)
#'
#' @param contig GRanges object with single contig
#' @param exons_gr GRanges object with reference exons
#' @param tolerance Maximum distance to consider as connected
#' @return data.table with matching exon or empty
find_upstream_connection <- function(contig, exons_gr, tolerance = 2) {
  
  chr <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  contig_start <- start(contig)
  contig_end <- end(contig)
  
  gene_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                         as.character(strand(exons_gr)) == contig_strand]
  
  if (contig_strand == "+") {
    nearby <- gene_exons[end(gene_exons) < contig_start &
                         end(gene_exons) >= (contig_start - tolerance)]
  } else {
    # biological 5' anchor on - strand is contig_end
    nearby <- gene_exons[start(gene_exons) > contig_end &
                         start(gene_exons) <= (contig_end + tolerance)]
  }
  
  if (length(nearby) == 0) return(data.table())
  
  closest_idx <- if (contig_strand == "+") which.max(end(nearby)) else which.min(start(nearby))
  closest <- nearby[closest_idx]
  
  result <- data.table(
    gene_id = closest$gene_id,
    transcript_id = closest$transcript_id,
    exon_boundary = if (contig_strand == "+") end(closest) else start(closest),
    distance = if (contig_strand == "+") {
      contig_start - end(closest)
    } else {
      start(closest) - contig_end
    }
  )
  
  return(result)
}

#' Find what exon the contig connects to at its 3' end (downstream)
#'
#' @param contig GRanges object with single contig
#' @param exons_gr GRanges object with reference exons
#' @param tolerance Maximum distance to consider as connected
#' @return data.table with matching exon or empty
find_downstream_connection <- function(contig, exons_gr, tolerance = 2) {
  
  chr <- as.character(seqnames(contig))
  contig_strand <- as.character(strand(contig))
  contig_start <- start(contig)
  contig_end <- end(contig)
  
  gene_exons <- exons_gr[as.character(seqnames(exons_gr)) == chr &
                         as.character(strand(exons_gr)) == contig_strand]
  
  if (contig_strand == "+") {
    nearby <- gene_exons[start(gene_exons) > contig_end &
                         start(gene_exons) <= (contig_end + tolerance)]
  } else {
    # biological 3' anchor on - strand is contig_start
    nearby <- gene_exons[end(gene_exons) < contig_start &
                         end(gene_exons) >= (contig_start - tolerance)]
  }
  
  if (length(nearby) == 0) return(data.table())
  
  closest_idx <- if (contig_strand == "+") which.min(start(nearby)) else which.max(end(nearby))
  closest <- nearby[closest_idx]
  
  result <- data.table(
    gene_id = closest$gene_id,
    transcript_id = closest$transcript_id,
    exon_boundary = if (contig_strand == "+") start(closest) else end(closest),
    distance = if (contig_strand == "+") {
      start(closest) - contig_end
    } else {
      contig_start - end(closest)
    }
  )
  
  return(result)
}

#' Check if a position matches an exon 5' boundary
#' For positive strand: exon start
#' For negative strand: exon end
is_exon_5_prime_boundary <- function(pos, exon, strand) {
  if (strand == "+") {
    return(pos == start(exon))
  } else {
    return(pos == end(exon))
  }
}

#' Check if a position matches an exon 3' boundary
#' For positive strand: exon end
#' For negative strand: exon start
is_exon_3_prime_boundary <- function(pos, exon, strand) {
  if (strand == "+") {
    return(pos == end(exon))
  } else {
    return(pos == start(exon))
  }
}

#' Classify a single exon contig (from exon_contigs_up)
#'
#' @param contig_idx Index of contig
#' @param contig_gr GRanges object with contig
#' @param exons_gr GRanges object with reference exons
#' @param contig_type "up" or "down"
#' @param tolerance Tolerance for coordinate matching
#' @return data.frame with classification
classify_exon_contig <- function(contig_idx, contig_gr, exons_gr,
                                 contig_type = "up", tolerance = 2) {
  
  contig <- contig_gr[contig_idx]
  
  result <- data.frame(
    contig_id = sprintf("%s:%d-%d:%s_%d",
                        seqnames(contig), start(contig), end(contig),
                        strand(contig), contig_idx),
    chr = as.character(seqnames(contig)),
    start = start(contig),
    end = end(contig),
    strand = as.character(strand(contig)),
    contig_type = contig_type,
    event_type = "unclassified",
    gene_id = NA_character_,
    transcript_id = NA_character_,
    confidence = NA_character_,
    stringsAsFactors = FALSE
  )
  
  region_info <- determine_region_type(contig, exons_gr)
  is_exonic <- region_info$is_exonic
  overlapping_exons <- region_info$overlapping_exons
  
  if (contig_type == "up") {
    
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
      # extends from 5' end of a reference exon → alternative_donor_late
      result$event_type  <- "alternative_donor_late"
      result$gene_id     <- upstream_conn$gene_id[1]
      result$transcript_id <- upstream_conn$transcript_id[1]
      result$confidence  <- "1"
    } else if (!has_upstream && has_downstream) {
      # extends from 3' end of a reference exon → alternative_acceptor_early
      result$event_type  <- "alternative_acceptor_early"
      result$gene_id     <- downstream_conn$gene_id[1]
      result$transcript_id <- downstream_conn$transcript_id[1]
      result$confidence  <- "1"
    } else if (has_upstream && has_downstream) {
      # connects to exons on both sides → intron retention
      result$event_type <- "intron_retention"
      result$gene_id    <- upstream_conn$gene_id[1]
      result$confidence <- "1"
    } else {
      result$event_type <- "cryptic_exon"
      result$confidence <- "2"
    }
    
  } else if (contig_type == "down") {
    
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
      result$event_type    <- "alternative_acceptor_late"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "1"
    } else if (!five_prime_match && three_prime_match) {
      result$event_type    <- "alternative_donor_early"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "1"
    } else if (!five_prime_match && !three_prime_match) {
      result$event_type    <- "exonic_cryptic"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "2"
    } else {
      result$event_type    <- "1_ES"
      result$gene_id       <- ref_exon$gene_id
      result$transcript_id <- ref_exon$transcript_id
      result$confidence    <- "3"
    }
  }
  
  return(result)
}

# ===== MAIN EXECUTION =====

cat("\n=== Loading Data ===\n")

# Load data
exons_gr <- load_reference(gtf_file)

cat("\n")
exon_contigs_up <- load_exon_contigs(exon_contigs_up_file, "up")
exon_contigs_down <- load_exon_contigs(exon_contigs_down_file, "down")

# Combine both sets
all_contigs_up <- exon_contigs_up
all_contigs_down <- exon_contigs_down

# ===== CLASSIFY CONTIGS =====

cat("\n=== Classifying Exon Contigs ===\n")

results <- list()
result_idx <- 1

# Classify up contigs
if (length(all_contigs_up) > 0) {
  cat(sprintf("\nClassifying %d up contigs...\n", length(all_contigs_up)))
  pb <- txtProgressBar(min = 0, max = length(all_contigs_up), style = 3)
  
  for (i in seq_along(all_contigs_up)) {
    results[[result_idx]] <- classify_exon_contig(i, all_contigs_up, exons_gr, 
                                                   "up", tolerance)
    result_idx <- result_idx + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

# Classify down contigs
if (length(all_contigs_down) > 0) {
  cat(sprintf("\nClassifying %d down contigs...\n", length(all_contigs_down)))
  pb <- txtProgressBar(min = 0, max = length(all_contigs_down), style = 3)
  
  for (i in seq_along(all_contigs_down)) {
    results[[result_idx]] <- classify_exon_contig(i, all_contigs_down, exons_gr,
                                                   "down", tolerance)
    result_idx <- result_idx + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

# ===== COMBINE RESULTS =====

cat("\n\n=== Combining Results ===\n")

if (length(results) == 0) {
  cat("No contigs to classify!\n")
  results_df <- data.frame()
} else {
  results_df <- bind_rows(results)
}

# ===== SUMMARY STATISTICS =====

cat("\n=== Classification Summary ===\n")

if (nrow(results_df) > 0) {
  print(table(results_df$event_type))
  
  cat("\nBreakdown by contig type:\n")
  print(table(results_df$contig_type, results_df$event_type))
} else {
  cat("No results to summarize\n")
}

# ===== WRITE OUTPUT =====

cat(sprintf("\nWriting results to %s\n", output_file))
write.table(
  results_df,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("\n✓ Exon contig classification complete!\n")

# Close log file
if (!is.null(log_file)) {
  sink()
  sink(type = "message")
  close(log_con)
}
