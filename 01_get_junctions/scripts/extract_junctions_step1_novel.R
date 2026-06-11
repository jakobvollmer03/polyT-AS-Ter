#!/usr/bin/env Rscript
## This script extracts all junctions (defined as gaps between exons) from a gtf of novel transcripts 
## (e. g. assembled by stringtie or manually curated). Snakemake will run extract_junctions_step1 to 
## get junctions from a reference gtf.


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
})

log_file <- snakemake@log[[1]]

# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}

# STEP 1: Extract Junctions from GTF Files
# This script extracts junctions and introns from both StringTie and reference GTFs
# and saves them as RDS files for later analysis

stringtie_gtf <- snakemake@input[["sample_gtf"]]


# Validate input files
if (!file.exists(stringtie_gtf)) {
  cat(sprintf("Error: StringTie GTF not found: %s\n", stringtie_gtf))
  quit(status = 1)
}


cat("\n=== STEP 1: Extract Novel Junctions ===\n")
cat(sprintf("StringTie GTF: %s\n", stringtie_gtf))

# ===== FUNCTION DEFINITIONS =====

#' Extract junctions from a GTF file
extract_junctions <- function(gtf) {
  
  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]
  
  if (length(exons) == 0) {
    warning("No exons found in GTF")
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
    warning("No novel junctions found")
    return(GRanges())
  }
  
  junctions <- do.call(c, junctions_list)
  
  # Remove duplicates
  junctions <- unique(junctions)
  
  return(junctions)
}

#' Extract introns from GTF
extract_introns <- function(gtf) {
  
  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]
  
  if (length(exons) == 0) {
    warning("No exons found in novel GTF")
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
    warning("No novel introns found")
    return(GRanges())
  }
  
  introns <- do.call(c, introns_list)
  
  # Remove duplicates
  introns <- unique(introns)
  
  return(introns)
}

#' Extract exons from GTF
extract_exons <- function(gtf) {
  
  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]
  
  if (length(exons) == 0) {
    warning("No exons found in novel GTF")
    return(GRanges())
  }
  
  # Remove duplicates
  exons <- unique(exons)
  
  return(exons)
}

# ===== MAIN PROCESSING =====

cat("\n[1/3] Loading StringTie GTF...\n")
stringtie <- import(stringtie_gtf)
cat(sprintf("      Loaded %d features\n", length(stringtie)))

cat("\n[2/3] Extracting junctions from StringTie...\n")
stringtie_junctions <- extract_junctions(stringtie)
cat(sprintf("      Found %d unique junctions\n", length(stringtie_junctions)))

cat("\n[3/3] Extracting exons from StringTie...\n")
stringtie_exons <- extract_exons(stringtie)
cat(sprintf("      Found %d unique exons\n", length(stringtie_exons)))

# ===== SAVE RESULTS =====

stringtie_jxn = snakemake@output[[1]]
stringtie_exn = snakemake@output[[2]]


cat("\n=== Saving Novel Ranges ===\n")
cat(sprintf("Saving StringTie junctions to: %s\n", stringtie_jxn))
saveRDS(stringtie_junctions, file = stringtie_jxn)

cat(sprintf("Saving StringTie exons to: %s\n", stringtie_exn))
saveRDS(stringtie_exons, file = stringtie_exn)

cat("\n=== Summary (Novel) ===\n")
cat(sprintf("StringTie junctions:  %d\n", length(stringtie_junctions)))
cat(sprintf("StringTie exons:      %d\n", length(stringtie_exons)))

# Close log file
if (!is.null(log_file)) {
  sink()
  sink(type = "message")
  close(log_con)
}