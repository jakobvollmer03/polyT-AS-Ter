#!/usr/bin/env Rscript
## This script extracts all junctions from a reference gtf for later comparison with novel transcripts for the
## purpose of identifying the exact novel splice junctions within the novel transcripts.


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# STEP 1: Extract Junctions from GTF Files
# and saves them as RDS files for later analysis

reference_gtf <- snakemake@input[["reference_gtf"]]
og_tx_ids <- snakemake@input[["og_id_file"]]
log_file <- snakemake@log[[1]]
# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}
# Validate input files
if (!file.exists(reference_gtf)) {
  cat(sprintf("Error: Reference GTF not found: %s\n", reference_gtf))
  quit(status = 1)
}

cat("\n=== STEP 1: Extract Reference Junctions ===\n")
cat(sprintf("Reference GTF: %s\n", reference_gtf))
cat(sprintf("Original Transcript IDs: %s\n", og_tx_ids))

# ===== FUNCTION DEFINITIONS =====

#' Extract junctions from a GTF file
extract_junctions <- function(gtf, og_tx_ids) {
  
  # Read original transcript IDs
  if (!is.null(og_tx_ids) && file.exists(og_tx_ids)) {
    og_tx_ids <- read_tsv(og_tx_ids, col_names = FALSE)$X1
  } else {
    og_tx_ids <- NULL
  }
  

  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]

  if (!is.null(og_tx_ids)) {
    exons <- exons[exons$transcript_id %in% og_tx_ids]
  }

  if (length(exons) == 0) {
    warning("No exons found in reference GTF")
    return(GRanges())
  }
  
  # Group exons by transcript
  exons_by_tx <- split(exons, exons$transcript_id)
  
  # Extract junctions from each transcript
  junctions_list <- lapply(exons_by_tx, function(tx_exons) {
    
    # Sort exons by start position
    tx_exons <- tx_exons[order(start(tx_exons))]
    
    # Need at least 2 exons to have junctions
    if (length(tx_exons) < 2) {
      return(NULL)
    }
    
    # Extract junction coordinates
    n_exons <- length(tx_exons)
    
    junction_starts <- end(tx_exons)[1:(n_exons-1)] + 1
    junction_ends <- start(tx_exons)[2:n_exons] - 1
    
    # Create GRanges for junctions
    junctions <- GRanges(
      seqnames = seqnames(tx_exons)[1],
      ranges = IRanges(start = junction_starts, end = junction_ends),
      strand = strand(tx_exons)[1]
    )
    
    return(junctions)
  })
  
  # Remove NULL entries
  junctions_list <- junctions_list[!sapply(junctions_list, is.null)]
  
  # Combine all junctions
  if (length(junctions_list) == 0) {
    warning("No reference junctions found")
    return(GRanges())
  }
  
  junctions <- do.call(c, junctions_list)
  
  # Remove duplicates
  junctions <- unique(junctions)
  
  return(junctions)
}

#' Extract introns from GTF
extract_introns <- function(gtf, og_tx_ids) {
  
  if (!is.null(og_tx_ids) && file.exists(og_tx_ids)) {
    og_tx_ids <- read_tsv(og_tx_ids, col_names = FALSE)$X1
  } else {
    og_tx_ids <- NULL
  }
  

  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]

  if (!is.null(og_tx_ids)) {
    exons <- exons[exons$transcript_id %in% og_tx_ids]
  }
  
  if (length(exons) == 0) {
    warning("No original exons found in reference GTF")
    return(GRanges())
  }
  
  # Group exons by transcript
  exons_by_tx <- split(exons, exons$transcript_id)
  
  # Extract introns from each transcript
  introns_list <- lapply(exons_by_tx, function(tx_exons) {
    
    # Sort exons by start position
    tx_exons <- tx_exons[order(start(tx_exons))]
    
    # Need at least 2 exons to have introns
    if (length(tx_exons) < 2) {
      return(NULL)
    }
    
    # Intron coordinates
    n_exons <- length(tx_exons)
    
    intron_starts <- end(tx_exons)[1:(n_exons-1)] + 1
    intron_ends <- start(tx_exons)[2:n_exons] - 1
    
    # Create GRanges for introns
    introns <- GRanges(
      seqnames = seqnames(tx_exons)[1],
      ranges = IRanges(start = intron_starts, end = intron_ends),
      strand = strand(tx_exons)[1]
    )
    
    return(introns)
  })
  
  # Remove NULL entries
  introns_list <- introns_list[!sapply(introns_list, is.null)]
  
  # Combine all introns
  if (length(introns_list) == 0) {
    warning("No reference introns found")
    return(GRanges())
  }
  
  introns <- do.call(c, introns_list)
  
  # Remove duplicates
  introns <- unique(introns)
  
  return(introns)
}

#' Extract exons from GTF
extract_exons <- function(gtf, og_tx_ids) {
  
  # Read original transcript IDs
  if (!is.null(og_tx_ids) && file.exists(og_tx_ids)) {
    og_tx_ids <- read_tsv(og_tx_ids, col_names = FALSE)$X1
  } else {
    og_tx_ids <- NULL
  }
  
  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]

  if (!is.null(og_tx_ids)) {
    exons <- exons[exons$transcript_id %in% og_tx_ids]
  }

  if (length(exons) == 0) {
    warning("No original exons found in reference GTF")
    return(GRanges())
  }
  
  # Remove duplicates
  exons <- unique(exons)
  
  return(exons)
}

# ===== MAIN PROCESSING =====


cat("\n[1/4] Loading reference GTF...\n")
reference <- import(reference_gtf)
cat(sprintf("     [1/4] Loaded %d features\n", length(reference)))

#if ("UCSC" %in% seqlevelsStyle(reference_gtf) ||
#    any(grepl("^chr", seqlevels(reference_gtf)))) {
#  
#  seqlevelsStyle(reference_gtf) <- "NCBI"
#  message("Converted chromosome names to NCBI style")
#}


cat("\n[2/4] Extracting junctions from reference...\n")
reference_junctions <- extract_junctions(reference, og_tx_ids)
cat(sprintf("     [2/4] Found %d unique junctions\n", length(reference_junctions)))

cat("\n[3/4] Extracting original introns from reference...\n")
reference_introns <- extract_introns(reference, og_tx_ids)
cat(sprintf("     [3/4] Found %d unique introns\n", length(reference_introns)))

cat("\n[4/4] Extracting original exons from reference...\n")
reference_exons <- extract_exons(reference, og_tx_ids)
cat(sprintf("     [4/4] Found %d unique exons\n", length(reference_exons)))

# ===== SAVE RESULTS =====

output_files <- c(
  reference_jxn = snakemake@output[[1]],
  reference_int = snakemake@output[[2]],
  reference_exn = snakemake@output[[3]]
)

cat("\n=== Saving Reference Ranges ===\n")

cat(sprintf("Saving reference junctions to: %s\n", output_files["reference_jxn"]))
saveRDS(reference_junctions, file = output_files["reference_jxn"])

cat(sprintf("Saving reference introns to:   %s\n", output_files["reference_int"]))
saveRDS(reference_introns, file = output_files["reference_int"])

cat(sprintf("Saving reference exons to:     %s\n", output_files["reference_exn"]))
saveRDS(reference_exons, file = output_files["reference_exn"])

cat("\n=== Summary ===\n")
cat(sprintf("Reference junctions:  %d\n", length(reference_junctions)))
cat(sprintf("Reference introns:    %d\n", length(reference_introns)))
cat(sprintf("Reference exons:      %d\n", length(reference_exons)))
