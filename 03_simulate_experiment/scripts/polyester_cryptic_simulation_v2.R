#!/usr/bin/env Rscript

# Polyester RNA-seq simulation for cryptic splice variants
# This script simulates RNA-seq data with:
# - Control: original transcripts expressed, novel transcripts barely/not expressed
# - Test: novel transcripts expressed, original transcripts barely/not expressed

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(polyester)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

# Get input/output from Snakemake
gtf_file <- snakemake@input[[1]]              # GTF file
fasta_file <- snakemake@input[[2]]            # Transcript sequences
og_transcripts_file <- snakemake@input[[3]]   # Original transcripts list
novel_transcripts_file <- snakemake@input[[4]] # Novel transcripts list

prefix <- snakemake@params[["prefix"]]
cfg <- snakemake@params[["cfg"]]
output_dir <- snakemake@params[["output_dir"]]
log_file <- snakemake@log[[1]]

# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}

# Read original and novel transcripts files
og_transcripts_file <- read.table(og_transcripts_file, header = TRUE)
novel_transcripts_file <- read.table(novel_transcripts_file, header = TRUE)

# Simulation parameters
num_reps <- cfg[["num_reps"]]           # Number of biological replicates per condition
reads_per_transcript <- cfg[["reads_per_transcript"]]  # Base number of reads per transcript
read_length <- cfg[["read_length"]]      # Read length (per mate for paired-end)
error_rate <- cfg[["error_rate"]]     # Sequencing error rate

# Expression fold changes
# For novel transcripts in test vs control
novel_fold_change_test <- cfg[["novel_fold_change_test"]]
novel_fold_change_control <- cfg[["novel_fold_change_control"]]     
# For original transcripts in test vs control  
original_fold_change_test <- cfg[["original_fold_change_test"]]
original_fold_change_control <- cfg[["original_fold_change_control"]] 

fraction_expressed <- cfg[["fraction_expressed"]] # fraction of regular genes expressed in the simulation (to avoid simulating reads for all transcripts, which is not realistic)

# ============================================================================
# FUNCTIONS
# ============================================================================

identify_novel_transcripts <- function(gtf, novel_ids, og_transcripts_ids) {
  # Extract exon entries and get unique transcript IDs
  # (GTF may only contain exon entries, not transcript entries)
  exons <- gtf[gtf$type == "exon"]
  
  if (length(exons) == 0) {
    # Fallback: if no exons, try to use transcript entries
    transcripts <- gtf[gtf$type == "transcript"]
    all_tx_ids <- transcripts$transcript_id
  } else {
    # Get unique transcript IDs from exon entries
    all_tx_ids <- unique(exons$transcript_id)
  }
  
  # Use provided novel IDs and original IDs to classify transcripts
  final_novel_ids <- all_tx_ids[all_tx_ids %in% novel_ids]
  final_original_ids <- all_tx_ids[all_tx_ids %in% og_transcripts_ids]
  final_regular_ids <- all_tx_ids[!all_tx_ids %in% novel_ids & !all_tx_ids %in% og_transcripts_ids]
   

  cat(sprintf("Found %d novel transcripts\n", length(final_novel_ids)))
  cat(sprintf("Found %d original transcripts\n", length(final_original_ids)))
  cat(sprintf("Found %d regular transcripts\n", length(final_regular_ids)))
  
  return(list(novel = final_novel_ids, original = final_original_ids, regular = final_regular_ids))
}

setup_fold_changes <- function(fasta_seqs, novel_ids, original_ids, regular_ids,
                                novel_fc_test, novel_fc_control, original_fc_test, original_fc_control,
                                reads_per_tx, num_reps = 1, fraction_expressed = 0.6, normalize_length = TRUE) {
  
  all_ids <- names(fasta_seqs)
  n_transcripts <- length(all_ids)
  
  # ===== LENGTH NORMALIZATION SETUP =====
  # Extract transcript lengths from FASTA sequences
  tx_lengths <- width(fasta_seqs)
  names(tx_lengths) <- all_ids
  
  # Use median as reference length (robust to outliers)
  median_length <- median(tx_lengths)
  
  cat("\n=== Transcript Length Statistics ===\n")
  cat(sprintf("Total transcripts: %d\n", n_transcripts))
  cat(sprintf("Length range: %d - %d bp\n", min(tx_lengths), max(tx_lengths)))
  cat(sprintf("Median length: %d bp (used as reference)\n", median_length))
  cat(sprintf("Mean length: %.0f bp\n", mean(tx_lengths)))
  
  if (normalize_length == TRUE) {
    # Calculate length normalization factors
    # Factor = transcript_length / median_length
    # Longer transcripts get more reads to achieve equal TPM
    length_factors <- tx_lengths / median_length
    cat("\nLength normalization: ENABLED\n")
    cat("All transcripts will have similar TPM values (not just read counts)\n")
  } else {
    # No normalization - all get same read count
    length_factors <- rep(1, n_transcripts)
    names(length_factors) <- all_ids
    cat("\nLength normalization: DISABLED\n")
    cat("All transcripts will get same read counts (TPM will vary by length)\n")
  }
  
  # ===== CREATE READS MATRIX =====
  reads_matrix <- matrix(0, nrow = n_transcripts, ncol = 2 * num_reps)
  rownames(reads_matrix) <- all_ids
  
  # Randomly select fraction of regular transcripts to be expressed
  n_regular <- length(regular_ids)
  n_expressed <- round(fraction_expressed * n_regular)
  expressed_regular <- sample(regular_ids, n_expressed)
  
  cat(sprintf("\nExpressed regular transcripts: %d / %d (%.1f%%)\n", 
              n_expressed, n_regular, 100 * fraction_expressed))
  
  # ===== ASSIGN READ COUNTS =====
  for (i in 1:n_transcripts) {
    tx_id <- all_ids[i]
    
    if (tx_id %in% novel_ids) {
      # Novel transcripts: low in control, high in test
      reads_matrix[i, 1:num_reps] <- novel_fc_control * reads_per_tx * length_factors[i]
      reads_matrix[i, (num_reps+1):(2*num_reps)] <- novel_fc_test * reads_per_tx * length_factors[i]
      
    } else if (tx_id %in% original_ids) {
      # Original transcripts: high in control, low in test
      reads_matrix[i, 1:num_reps] <- original_fc_control * reads_per_tx * length_factors[i]
      reads_matrix[i, (num_reps+1):(2*num_reps)] <- original_fc_test * reads_per_tx * length_factors[i]
      
    } else if (tx_id %in% expressed_regular) {
      # Regular (expressed): uniform across both conditions
      reads_matrix[i, 1:num_reps] <- reads_per_tx * length_factors[i]
      reads_matrix[i, (num_reps+1):(2*num_reps)] <- reads_per_tx * length_factors[i]
      
    }
    # else: not expressed, remains 0
  }
  
  # ===== DIAGNOSTIC OUTPUT =====
  cat("\n=== Read Count Examples ===\n")
  
  # Show examples from each category
  if (length(expressed_regular) > 0) {
    sample_idx <- head(which(all_ids %in% expressed_regular), 3)
    cat("\nRegular transcripts (should have uniform TPM):\n")
    for (idx in sample_idx) {
      cat(sprintf("  %s: length=%d bp, factor=%.2f, control_reads=%.0f, test_reads=%.0f\n",
                  all_ids[idx], tx_lengths[idx], length_factors[idx],
                  reads_matrix[idx, 1], reads_matrix[idx, num_reps+1]))
    }
  }
  
  if (any(all_ids %in% novel_ids)) {
    sample_idx <- head(which(all_ids %in% novel_ids), 2)
    cat("\nNovel transcripts (low → high):\n")
    for (idx in sample_idx) {
      cat(sprintf("  %s: length=%d bp, control_reads=%.0f, test_reads=%.0f (%.0fx change)\n",
                  all_ids[idx], tx_lengths[idx],
                  reads_matrix[idx, 1], reads_matrix[idx, num_reps+1],
                  reads_matrix[idx, num_reps+1] / max(reads_matrix[idx, 1], 1)))
    }
  }
  
  if (any(all_ids %in% original_ids)) {
    sample_idx <- head(which(all_ids %in% original_ids), 2)
    cat("\nOriginal transcripts (high → low):\n")
    for (idx in sample_idx) {
      cat(sprintf("  %s: length=%d bp, control_reads=%.0f, test_reads=%.0f (%.0fx change)\n",
                  all_ids[idx], tx_lengths[idx],
                  reads_matrix[idx, 1], reads_matrix[idx, num_reps+1],
                  reads_matrix[idx, num_reps+1] / max(reads_matrix[idx, 1], 1)))
    }
  }
  
  cat("\n")
  
  return(list(
    reads_matrix = reads_matrix, 
    expressed_regular = expressed_regular,
    tx_lengths = tx_lengths,
    median_length = median_length
  ))
}
# ============================================================================
# MAIN EXECUTION
# ============================================================================

tryCatch({
  cat("=== Polyester Cryptic Splice Simulation ===\n\n")
  
  # Load GTF file
  cat("Loading GTF file...\n")
  gtf <- import(gtf_file)
  
  # Identify novel vs original transcripts
  cat("Identifying novel and original transcripts...\n")
  transcript_groups <- identify_novel_transcripts(gtf, novel_transcripts_file$transcript_id, 
                                                   og_transcripts_file$transcript_id)
  
  # Load transcript sequences
  cat("Loading transcript sequences...\n")
  fasta_seqs <- readDNAStringSet(fasta_file)
  names(fasta_seqs) <- sapply(strsplit(names(fasta_seqs), " "), `[`, 1)
  
  cat(sprintf("Loaded %d transcript sequences\n\n", length(fasta_seqs)))
  
  cat("Sample FASTA IDs:\n")
  print(head(names(fasta_seqs)))
  cat("\nSample novel IDs from GTF:\n")
  print(head(transcript_groups$novel))
  cat("\nSample original IDs from GTF:\n")
  print(head(transcript_groups$original))
  cat("\nDo any FASTA IDs match novel IDs?", any(names(fasta_seqs) %in% transcript_groups$novel), "\n")
  cat("Do any FASTA IDs match original IDs?", any(names(fasta_seqs) %in% transcript_groups$original), "\n\n")
  cat(sprintf("Loaded %d transcript sequences\n\n", length(fasta_seqs)))
  
  # Set up fold changes
  cat("Setting up fold changes...\n")
  fc_setup <- setup_fold_changes(
    fasta_seqs, 
    transcript_groups$novel, 
    transcript_groups$original,
    transcript_groups$regular,
    novel_fold_change_test,
    novel_fold_change_control,
    original_fold_change_test,
    original_fold_change_control,
    reads_per_transcript,
    num_reps  
  )
  
  reads_matrix <- fc_setup$reads_matrix
  
  # Validation check
  novel_idx <- which(rownames(reads_matrix) %in% transcript_groups$novel)[1]
  original_idx <- which(rownames(reads_matrix) %in% transcript_groups$original)[1]
  # report paths of input files for debugging
  cat("Input files:\n")
  cat(sprintf("  GTF file: %s\n", gtf_file))
  cat(sprintf("  FASTA file: %s\n", fasta_file))
  cat(sprintf("  Original transcripts file: %s\n", og_transcripts_file))
  cat(sprintf("  Novel transcripts file: %s\n", novel_transcripts_file))
  
  cat(sprintf("Sample novel transcript: Control=%.2f, Test=%.2f\n", 
              reads_matrix[novel_idx, 1], reads_matrix[novel_idx, num_reps + 1]))
  cat(sprintf("Sample original transcript: Control=%.2f, Test=%.2f\n", 
              reads_matrix[original_idx, 1], reads_matrix[original_idx, num_reps + 1]))
  
  cat(sprintf("Regular transcripts expressed: %d (%.1f%%)\n", 
              length(fc_setup$expressed_regular),
              100 * length(fc_setup$expressed_regular) / length(transcript_groups$regular)))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run Polyester simulation
  cat("Running Polyester simulation...\n")
  cat(sprintf("  Replicates per condition: %d\n", num_reps))
  cat(sprintf("  Read length: %d bp\n", read_length))
  cat(sprintf("  Reads per transcript: ~%d\n", reads_per_transcript))
  cat(sprintf("  Error rate: %f\n\n", error_rate))
  
  simulate_experiment_countmat(
    fasta = fasta_file,
    outdir = output_dir,
    num_reps = c(num_reps, num_reps),
    readmat = reads_matrix,
    readlen = read_length,
    error_rate = error_rate,
    paired = TRUE,
    strand_specific = TRUE
  )
  
  # Rename Polyester output files to include prefix
  cat("Renaming FASTA files with prefix...\n")
  for (i in 1:(2*num_reps)) {
    sample_num <- sprintf("%02d", i)
    old_name_1 <- file.path(output_dir, paste0("sample_", sample_num, "_1.fasta"))
    old_name_2 <- file.path(output_dir, paste0("sample_", sample_num, "_2.fasta"))
    new_name_1 <- file.path(output_dir, paste0(prefix, "_sample_", sample_num, "_1.fasta"))
    new_name_2 <- file.path(output_dir, paste0(prefix, "_sample_", sample_num, "_2.fasta"))
    
    if (file.exists(old_name_1)) {
      file.rename(old_name_1, new_name_1)
    }
    if (file.exists(old_name_2)) {
      file.rename(old_name_2, new_name_2)
    }
  }
  
  # Save simulation metadata
  cat("Saving simulation metadata...\n")
  metadata <- data.frame(
    sample_id = paste0("sample_", sprintf("%02d", 1:(2*num_reps))),
    condition = rep(c("control", "test"), each = num_reps),
    replicate = rep(1:num_reps, 2),
    fastq_1 = file.path(output_dir, paste0(prefix, "_sample_", sprintf("%02d", 1:(2*num_reps)), "_1.fasta")),
    fastq_2 = file.path(output_dir, paste0(prefix, "_sample_", sprintf("%02d", 1:(2*num_reps)), "_2.fasta"))
  )
  
  write.csv(metadata, file.path(output_dir, paste0(prefix, "_sample_metadata.csv")), row.names = FALSE)
  
  # Save transcript group assignments
  transcript_assignments <- data.frame(
    transcript_id = c(transcript_groups$novel, transcript_groups$original),
    group = c(rep("novel", length(transcript_groups$novel)),
              rep("original", length(transcript_groups$original)))
  )
  write.csv(transcript_assignments, 
            file.path(output_dir, paste0(prefix, "_transcript_groups.csv")), 
            row.names = FALSE)
  
  # Save fold change matrix for reference
  write.csv(reads_matrix, 
            file.path(output_dir, paste0(prefix, "_reads_matrix.csv")))
  
  cat("\n=== Simulation complete! ===\n")
  cat(sprintf("Output directory: %s\n", output_dir))
  cat(sprintf("  - FASTA files: %s_sample_*.fasta\n", prefix))
  cat(sprintf("  - Metadata: %s_sample_metadata.csv\n", prefix))
  cat(sprintf("  - Transcript groups: %s_transcript_groups.csv\n", prefix))
  cat(sprintf("  - Read Matrix: %s_reads_matrix.csv\n", prefix))
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("Backtrace:\n")
  print(traceback())
  quit(status = 1)
})

# Close log file
if (!is.null(log_file)) {
  sink()
  sink(type = "message")
  close(log_con)
}