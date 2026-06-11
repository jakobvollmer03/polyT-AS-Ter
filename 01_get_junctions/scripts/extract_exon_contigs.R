suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(tools)
})
log_file <- snakemake@log[[1]]
# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}
input_files <- c(
  sample_exons = snakemake@input[[1]],
  reference_exons = snakemake@input[[2]],
  og_id_file = snakemake@input[[3]]
)
output_file_up <- snakemake@output[[1]]
output_file_down <- snakemake@output[[2]]
prefix <- snakemake@params[["prefix"]]

for (fname in input_files) {
  if (!file.exists(fname)) {
    cat(sprintf("Error: Required file not found: %s\n", fname))
    cat("Make sure you ran Step 1 with the same prefix.\n")
    quit(status = 1)
  }
}

cat("\n=== STEP 3: Extract Exon Contigs ===\n")
cat(sprintf("Input prefix:     %s\n", prefix))
cat(sprintf("Sample exons:     %s\n", input_files["sample_exons"]))
cat(sprintf("Reference exons:  %s\n", input_files["reference_exons"]))
cat(sprintf("Output file (up): %s\n", output_file_up))
cat(sprintf("Output file (down): %s\n\n", output_file_down))

# Find sample exons that are not identical to reference exons
cat("Loading exon data...\n")
sample_exons <- readRDS(input_files["sample_exons"])
reference_exons <- readRDS(input_files["reference_exons"])
og_id_df <- read_tsv(input_files["og_id_file"], col_types = cols())

cat(sprintf("Loaded %d sample exons\n", length(sample_exons)))
cat(sprintf("Loaded %d reference exons\n", length(reference_exons)))
cat(sprintf("Loaded %d original transcript IDs\n", nrow(og_id_df)))

# Filter reference exons to only include those from original transcripts
# (the transcripts that the novel transcripts are derived from)
og_exons <- reference_exons[reference_exons$transcript_id %in% og_id_df$transcript_id]
cat(sprintf("Filtered to %d exons from original transcripts\n", length(og_exons)))

# ===== FUNCTION DEFINITIONS =====

#' Extract exon segments not in a reference set
#' For each exon in query_exons, find parts that don't overlap with any exon in ref_exons
#' (considering same chromosome and strand)
extract_novel_segments <- function(query_exons, ref_exons, description = "") {
  
  if (length(query_exons) == 0) {
    cat(sprintf("No %s to process\n", description))
    return(GRanges())
  }
  
  # Standardize seqlevels to avoid factor level mismatches
  all_chroms <- unique(c(seqlevels(query_exons), seqlevels(ref_exons)))
  seqlevels(query_exons) <- all_chroms
  seqlevels(ref_exons) <- all_chroms
  
  novel_segments <- list()
  
  # For each query exon, find what parts don't overlap with reference
  for (i in seq_along(query_exons)) {
    query_ex <- query_exons[i]
    query_chr <- as.character(seqnames(query_ex))
    query_strand <- as.character(strand(query_ex))
    
    # Find reference exons on same chromosome and strand
    matching_ref <- ref_exons[as.character(seqnames(ref_exons)) == query_chr & 
                               as.character(strand(ref_exons)) == query_strand]
    
    if (length(matching_ref) == 0) {
      # No overlapping reference exons - entire query exon is novel
      novel_segments[[i]] <- query_ex
    } else {
      # Find overlaps
      hits <- findOverlaps(query_ex, matching_ref, ignore.strand = FALSE)
      
      if (length(hits) == 0) {
        # Query exon doesn't overlap any reference exon
        novel_segments[[i]] <- query_ex
      } else {
        # Get the union of overlapping reference exons
        overlapping_refs <- matching_ref[subjectHits(hits)]
        ref_union <- reduce(overlapping_refs)
        
        # Subtract overlapping regions from query exon
        complement <- setdiff(query_ex, ref_union)
        
        if (length(complement) > 0) {
          novel_segments[[i]] <- complement
        }
      }
    }
  }
  
  # Combine all novel segments
  if (length(novel_segments) == 0) {
    return(GRanges())
  }
  
  # Filter out NULL entries and combine
  novel_segments <- novel_segments[!sapply(novel_segments, is.null)]
  
  if (length(novel_segments) == 0) {
    return(GRanges())
  }
  
  combined <- do.call(c, novel_segments)
  
  # Standardize seqlevels for output
  combined_chroms <- unique(seqlevels(combined))
  seqlevels(combined) <- combined_chroms
  
  return(combined)
}

# ===== PROCESS EXONS =====

cat("\n=== Processing Sample vs Reference ===\n")

cat("\n[1/2] Finding sample exon segments not in reference...\n")
sample_novel <- extract_novel_segments(sample_exons, og_exons, 
                                       "sample exons")
sample_novel <- unique(sample_novel)
cat(sprintf("      Found %d unique sample segments\n", length(sample_novel)))

cat("\n[2/2] Finding reference exon segments not in sample...\n")
reference_novel <- extract_novel_segments(og_exons, sample_exons,
                                          "original exons")
reference_novel <- unique(reference_novel)
cat(sprintf("      Found %d unique reference segments\n", length(reference_novel)))

# ===== FORMAT OUTPUT =====

cat("\n=== Preparing Output ===\n")

# Convert to data frames
sample_novel_df <- data.frame(
  chr = as.character(seqnames(sample_novel)),
  start = start(sample_novel),
  end = end(sample_novel),
  strand = as.character(strand(sample_novel)),
  stringsAsFactors = FALSE
)

reference_novel_df <- data.frame(
  chr = as.character(seqnames(reference_novel)),
  start = start(reference_novel),
  end = end(reference_novel),
  strand = as.character(strand(reference_novel)),
  stringsAsFactors = FALSE
)

# Sort by chromosome and position
if (nrow(sample_novel_df) > 0) {
  sample_novel_df <- sample_novel_df %>% arrange(chr, start, end)
}

if (nrow(reference_novel_df) > 0) {
  reference_novel_df <- reference_novel_df %>% arrange(chr, start, end)
}

# ===== WRITE OUTPUT =====

cat("\n=== Writing Output Files ===\n")

cat(sprintf("Writing sample-specific exon segments (up) to: %s\n", output_file_up))
write.table(
  sample_novel_df,
  file = output_file_up,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(sprintf("Writing reference-specific exon segments (down) to: %s\n", output_file_down))
write.table(
  reference_novel_df,
  file = output_file_down,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ===== SUMMARY =====

cat("\n=== Summary ===\n")
cat(sprintf("Sample-specific exon segments (up):      %d\n", nrow(sample_novel_df)))
cat(sprintf("Reference-specific exon segments (down): %d\n", nrow(reference_novel_df)))

cat("\n✓ Step 3 complete!\n")
cat("\nOutput files created:\n")
cat(sprintf("  - %s (sample exons or parts not in reference)\n", output_file_up))
cat(sprintf("  - %s (reference exons or parts not in sample)\n", output_file_down))