#!/usr/bin/env Rscript

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
# Get paths from Snakemake params and input
stringtie_gtf <- snakemake@params[["sample_gtf"]]
output_file <- snakemake@output[[1]]
prefix <- snakemake@params[["prefix"]]
# Define expected input files (from Step 1 outputs)
input_files <- c(
  stringtie_jxn = snakemake@input[[1]],  # from Step 1
  reference_jxn = snakemake@input[[2]],  # from Step 1
  reference_int = snakemake@input[[3]],   # from Step 1
  reference_exons = snakemake@input[[4]]   # from Step 1
)

# Validate input files exist
for (fname in input_files) {
  if (!file.exists(fname)) {
    cat(sprintf("Error: Required file not found: %s\n", fname))
    cat("Make sure you ran Step 1 with the same prefix.\n")
    quit(status = 1)
  }
}

if (!file.exists(stringtie_gtf)) {
  cat(sprintf("Error: StringTie GTF not found: %s\n", stringtie_gtf))
  quit(status = 1)
}

cat("\n=== STEP 2: Identify Novel Junctions ===\n")
cat(sprintf("Input prefix:     %s\n", prefix))
cat(sprintf("StringTie GTF:    %s\n", stringtie_gtf))
cat(sprintf("Output file:      %s\n\n", output_file))

# Load libraries
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tools)
})

# ===== FUNCTION DEFINITIONS =====

#' Find retained introns
#' @param stringtie_gtf GRanges object with StringTie GTF
#' @param reference_introns GRanges object with reference introns
#' @param reference_exons GRanges object with reference exons
#' @return GRanges object with retained intron coordinates
find_retained_introns <- function(stringtie_gtf, reference_introns, reference_exons) {

  ## Normalize introns
  if (is(reference_introns, "GRangesList")) {
    reference_introns <- unlist(reference_introns, use.names = FALSE)
  }

  if (!is(reference_introns, "GRanges")) {
    stop("reference_introns must be GRanges or GRangesList")
  }

  if (!is(stringtie_gtf, "GRanges")) {
    stop("stringtie_gtf must be a GRanges")
  }

  if (length(reference_introns) == 0) {
    return(GRanges())
  }

  ## Extract StringTie exons
  stringtie_exons <- stringtie_gtf[stringtie_gtf$type == "exon"]

  if (length(stringtie_exons) == 0) {
    return(GRanges())
  }
  stringtie_exons_filtered <- stringtie_exons[!overlapsAny(stringtie_exons, reference_exons, type = "equal")]

  ## Find introns fully contained within exons
  hits <- findOverlaps(
    reference_introns,
    stringtie_exons_filtered,
    type = "within",
    ignore.strand = FALSE
  )

  if (length(hits) == 0) {
    return(GRanges())
  }

  ## Return unique retained introns
  retained_introns <- unique(
    reference_introns[queryHits(hits)]
  )

  return(retained_introns)
}


# ===== LOAD SAVED DATA =====

cat("\n[1/4] Loading StringTie junctions...\n")
stringtie_junctions <- readRDS(input_files["stringtie_jxn"])
cat(sprintf("      Loaded %d junctions\n", length(stringtie_junctions)))
#stringtie_junctions <- do.call(c, stringtie_junctions) # convert list to GRanges

cat("\n[2/4] Loading reference transcripts...\n")
reference_junctions <- readRDS(input_files["reference_jxn"])
cat(sprintf("      Loaded %d transcripts\n", length(reference_junctions)))
#reference_junctions <- do.call(c, reference_junctions) # convert list to GRanges

cat("\n[3/4] Loading reference introns...\n")
reference_introns <- readRDS(input_files["reference_int"])
cat(sprintf("      Loaded %d introns\n", length(reference_introns)))
#reference_introns <- do.call(c, reference_introns) # convert list to GRanges

cat("\n[4/4] Loading reference exons...\n")
reference_exons <- readRDS(input_files["reference_exons"])
cat(sprintf("      Loaded %d exons\n", length(reference_exons)))
#reference_exons <- do.call(c, reference_exons) # convert list to GRanges

cat("\n[5/5] Loading StringTie GTF (for retained intron detection)...\n")
stringtie <- import(stringtie_gtf)
cat(sprintf("      Loaded %d features\n", length(stringtie)))

# ===== IDENTIFY NOVEL JUNCTIONS =====

cat("\n=== Identifying Novel Features ===\n")

cat("\n[1/2] Identifying novel junctions...\n")
# Find junctions in StringTie that don't exist in reference
# Use overlapsAny with type="equal" for exact matching
stringtie_junctions <- do.call(c, stringtie_junctions) # convert list to GRanges
reference_junctions <- do.call(c, reference_junctions) # convert list to GRanges

novel_junctions <- stringtie_junctions[!overlapsAny(stringtie_junctions, 
                                                      reference_junctions, 
                                                      type = "equal")]
novel_junctions <- unique(novel_junctions)
cat(sprintf("      Found %d novel junctions\n", length(novel_junctions)))

reference_introns <- do.call(c, reference_introns) # convert list to GRanges
cat("\n[2/2] Identifying retained introns...\n")
retained_introns <- find_retained_introns(stringtie, reference_introns, reference_exons)
retained_introns <- unique(retained_introns)
cat(sprintf("      Found %d retained introns\n", length(retained_introns)))

# ===== DIAGNOSTIC: Check for consecutive/overlapping retained introns =====
if (length(retained_introns) > 1) {
  # Sort by chromosome and position
  retained_sorted <- retained_introns[order(seqnames(retained_introns), 
                                             start(retained_introns))]
  
  # Find overlaps or adjacency between retained introns
  overlaps <- findOverlaps(retained_sorted, retained_sorted, 
                          ignore.strand = FALSE)
  
  # Remove self-overlaps
  overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
  
  if (length(overlaps) > 0) {
    n_overlapping <- length(unique(queryHits(overlaps)))
    cat(sprintf("\n      NOTE: Found %d retained introns that overlap or are adjacent to other retained introns\n", 
                n_overlapping))
    cat("      Each intron is reported separately (not merged) for granular downstream analysis.\n")
    cat("      This may represent:\n")
    cat("        - Multiple independent retention events\n")
    cat("        - A large retention spanning multiple annotated introns\n")
    cat("        - Assembly or annotation artifacts\n")
  }
}

# ===== COMBINE AND FORMAT OUTPUT =====

cat("\n=== Preparing Output ===\n")

# Combine novel junctions and retained introns
# NOTE: c() on GRanges concatenates without merging
all_novel <- c(novel_junctions, retained_introns)

# Remove exact duplicates only - this will NOT merge adjacent ranges
all_novel <- unique(all_novel)

cat(sprintf("Combined features: %d (from %d junctions + %d introns)\n", 
            length(all_novel), length(novel_junctions), length(retained_introns)))

# Verify no unexpected merging occurred
expected_count <- length(unique(c(novel_junctions, retained_introns)))
if (length(all_novel) != expected_count) {
  warning(sprintf("Unexpected change in count: expected %d, got %d", 
                  expected_count, length(all_novel)))
}

# Convert to data frame


novel_junctions_df <- data.frame(
  chr = as.character(seqnames(novel_junctions)),
  start = start(novel_junctions),
  end = end(novel_junctions),
  strand = as.character(strand(novel_junctions)),
  type = "splice_junction",
  stringsAsFactors = FALSE
)

retained_introns_df <- data.frame(
  chr = as.character(seqnames(retained_introns)),
  start = start(retained_introns),
  end = end(retained_introns),
  strand = as.character(strand(retained_introns)),
  type = "intron_retention",
  stringsAsFactors = FALSE
)

result_df <- data.frame(
  chr = as.character(seqnames(all_novel)),
  start = start(all_novel),
  end = end(all_novel),
  strand = as.character(strand(all_novel)),
  type = ifelse(all_novel %in% novel_junctions, "splice_junction", "intron_retention"),
  stringsAsFactors = FALSE
)
# Sort by chromosome and position
result_df <- result_df %>% arrange(chr, start, end)
novel_junctions_df <- novel_junctions_df %>% arrange(chr, start, end)
retained_introns_df <- retained_introns_df %>% arrange(chr, start, end)

# ===== WRITE OUTPUT =====

cat("\n=== Writing Output Files ===\n")

output_file <- snakemake@output[[1]]
jxn_file <- snakemake@output[[2]]
int_file <- snakemake@output[[3]]

cat(sprintf("Writing combined novel features to: %s\n", output_file))
write.table(
  result_df,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(sprintf("Writing novel junctions to: %s\n", jxn_file))
write.table(
  novel_junctions_df,
  file = jxn_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(sprintf("Writing retained introns to: %s\n", int_file))
write.table(
  retained_introns_df,
  file = int_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ===== SUMMARY =====

cat("\n=== Summary ===\n")
cat(sprintf("Total novel features:    %d\n", nrow(result_df)))
cat(sprintf("  Novel junctions:       %d\n", nrow(novel_junctions_df)))
cat(sprintf("  Retained introns:      %d\n", nrow(retained_introns_df)))

cat("\n✓ Step 2 complete!\n")
cat("\nOutput files created:\n")
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", jxn_file))
cat(sprintf("  - %s\n", int_file))
