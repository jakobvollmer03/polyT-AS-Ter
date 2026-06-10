suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(here)
})
source(here::here("scripts/utils.R"))

#' Normalise FASTA sequence names to bare chromosome identifiers.
#'
#' Ensembl FASTA headers look like:
#'   "1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF"
#' The chromosome name is always the first whitespace-delimited token.
#' The GTF (and therefore the modifications table) already uses bare names.
#' This function strips everything from the first space onward so lookups work.
#'
#' @param fasta  DNAStringSet as returned by readDNAStringSet
#' @return DNAStringSet with names replaced by first-token chromosome identifiers
normalise_fasta_names <- function(fasta) {
  raw_names  <- names(fasta)
  norm_names <- sub("\\s.*$", "", raw_names)   # keep only up to first whitespace

  # Guard: normalisation must produce unique names.
  # If two headers share the same first token the FASTA is malformed.
  if (anyDuplicated(norm_names)) {
    dups <- norm_names[duplicated(norm_names)]
    stop(sprintf(
      "normalise_fasta_names: normalised chromosome names are not unique. ",
      "Duplicated names: %s. ",
      "This suggests the FASTA has malformed headers — each sequence must have ",
      "a unique identifier as its first whitespace-delimited token.",
      paste(unique(dups), collapse = ", ")))
  }

  if (!identical(raw_names, norm_names)) {
    message(sprintf(
      "[INFO] Normalised %d FASTA sequence names (stripped Ensembl header metadata). Example: '%s' → '%s'",
      sum(raw_names != norm_names),
      raw_names[1L],
      norm_names[1L]
    ))
  }

  names(fasta) <- norm_names
  fasta
}

#' Check that no modification position falls within an exon of a different transcript.
#'
#' Stops with an informative error if any mod position overlaps a reference exon
#' belonging to a transcript other than the one that generated the modification.
#'
#' @param mods   data.frame with columns chr, pos1, pos2, transcript_id
#' @param exons  GRanges of exon features from the collapsed reference
#check_mod_safety <- function(mods, exons) {
#  # Ensure only exon-type features are checked (belt-and-suspenders)
#  if (!is.null(exons$type)) exons <- exons[exons$type == "exon"]
#
#  exons_dt <- data.table(
#    chr           = as.character(seqnames(exons)),
#    start         = start(exons),
#    end           = end(exons),
#    transcript_id = exons$transcript_id
#  )
#
#  conflicts <- 0L
#  for (i in seq_len(nrow(mods))) {
#    m    <- mods[i, ]
#    hits <- exons_dt[
#      chr == m$chr &
#      ((start <= m$pos1 & end >= m$pos1) |
#       (start <= m$pos2 & end >= m$pos2)) &
#      transcript_id != m$transcript_id
#    ]
#    if (nrow(hits) > 0L) {
#      # For early-donor / late-acceptor events the modified positions are within
#      # the canonical exon of the source transcript — this is expected and safe
#      # (both the source and the novel transcript are grounded in the same
#      # modified FASTA). Only flag positions that are exonic in a *different* gene.
#      other_gene_hits <- hits[
#        !grepl(paste0("^", gsub("_novel.*$", "", m$transcript_id)), transcript_id)
#      ]
#      if (nrow(other_gene_hits) > 0L) {
#        conflicts <- conflicts + 1L
#        log_warn(sprintf(
#          "FASTA mod at %s:%d-%d (transcript %s) overlaps exon of unrelated transcript %s. ",
#          "This candidate should be excluded — adjust seed or tighten filters.",
#          m$chr, m$pos1, m$pos2, m$transcript_id, other_gene_hits$transcript_id[1L]))
#      }
#    }
#  }
#  if (conflicts > 0L)
#    stop(sprintf(
#      "check_mod_safety: %d modification(s) conflict with unrelated exons. ",
#      "See warnings above. Re-run with a different seed or stricter overlap filters.",
#      conflicts))
#  invisible(NULL)
#}

#' Apply a table of single-base substitutions to a DNAStringSet.
#'
#' Each row of mods specifies one pair of bases to write (pos1/base1, pos2/base2).
#' Positions are 1-based.  Operates via character conversion for safety.
#'
#' @param fasta  DNAStringSet (names must already be normalised to bare chr IDs)
#' @param mods   data.frame with columns chr, pos1, pos2, base1, base2
#' @return Modified DNAStringSet (same names as input)
apply_fasta_mods <- function(fasta, mods) {
  seqs <- as.character(fasta)   # named character vector; safe for substr<- assignment

  fasta_chroms   <- names(seqs)
  mod_chroms     <- unique(mods$chr)
  missing_chroms <- setdiff(mod_chroms, fasta_chroms)

  if (length(missing_chroms) > 0L) {
    stop(sprintf(
      "apply_fasta_mods: chromosome name mismatch after normalisation. ",
      "Modifications reference: %s. ",
      "FASTA contains: %s. ",
      "If the FASTA uses a non-Ensembl naming convention (e.g. 'chr1' vs '1'), ",
      "add a manual name mapping before calling apply_fasta_mods.",
      paste(missing_chroms, collapse = ", "),
      paste(head(fasta_chroms, 10L), collapse = ", ")))
  }

  for (i in seq_len(nrow(mods))) {
    m        <- mods[i, ]
    chr_name <- m$chr
    substr(seqs[chr_name], m$pos1, m$pos1) <- m$base1
    substr(seqs[chr_name], m$pos2, m$pos2) <- m$base2
  }
  DNAStringSet(seqs)
}

#' Convert 1-based modification table to BED format (0-based half-open intervals).
#'
#' Each row of mods becomes one BED record spanning pos1 to pos2 (inclusive in
#' 1-based coords → [pos1-1, pos2) in BED).
#'
#' @param mods  data.frame with columns chr, pos1, pos2, transcript_id, site
#' @return data.frame with BED columns: chrom, chromStart, chromEnd, name, score, strand
mods_to_bed <- function(mods) {
  data.frame(
    chrom      = mods$chr,
    chromStart = mods$pos1 - 1L,   # 1-based → 0-based
    chromEnd   = mods$pos2,        # half-open: pos2 (1-based) = pos2 (0-based end)
    name       = paste0(mods$transcript_id, "_", mods$site),
    score      = 0L,
    strand     = ".",
    stringsAsFactors = FALSE
  )
}

# ── Snakemake glue ────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")

  exons <- load_gtf(snakemake@input[["collapsed_gtf"]])
  if (!is.null(exons$type)) exons <- exons[exons$type == "exon"]

  mods <- read.table(snakemake@input[["mods"]], header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)

  #log_info("Checking modification safety against collapsed reference exons...")
  #check_mod_safety(mods, exons)

  log_info("Loading genome FASTA...")
  fasta <- readDNAStringSet(snakemake@input[["fasta"]])

  # Normalise Ensembl-style headers: "1 dna:chromosome ..." → "1"
  # This must happen before apply_fasta_mods so chromosome name lookups succeed.
  fasta <- normalise_fasta_names(fasta)

  log_info(sprintf("Applying %d base substitutions...", nrow(mods)))
  fasta <- apply_fasta_mods(fasta, mods)
  writeXStringSet(fasta, snakemake@output[["fasta"]])

  bed <- mods_to_bed(mods)
  write.table(bed, snakemake@output[["bed"]],
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  log_info(sprintf("Done. Modified %d positions across %d novel junctions.",
                   nrow(mods) * 2L, nrow(mods) / 2L))

  sink(); sink(type = "message")
  close(log_con)
}