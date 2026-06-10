suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(truncnorm)
})

#' Import a GTF and return exon features only
load_gtf <- function(path) {
  gtf <- rtracklayer::import(path)
  if (!"type" %in% names(mcols(gtf))) stop("GTF has no 'type' column: ", path)
  exons <- gtf[gtf$type == "exon"]
  if (length(exons) == 0L) stop("No exon features found in GTF: ", path)
  exons
}

#' Build exon boundary index for fast lookup (adapted from classifier — copy-before-setkey order corrected vs. original)
build_exon_index <- function(exons) {
  dt <- data.table(
    chr          = as.character(seqnames(exons)),
    strand       = as.character(strand(exons)),
    exon_start   = start(exons),
    exon_end     = end(exons),
    gene_id      = exons$gene_id,
    transcript_id = exons$transcript_id,
    exon_number  = exons$exon_number
  )
  end_index   <- copy(dt); setkey(end_index,   chr, strand, exon_end)
  start_index <- copy(dt); setkey(start_index, chr, strand, exon_start)
  list(exon_end = end_index, exon_start = start_index)
}

#' Total spliced length for one transcript
spliced_length <- function(exons_gr, tx_id) {
  hits <- exons_gr[exons_gr$transcript_id == tx_id]
  if (length(hits) == 0L) stop("transcript_id not found: ", tx_id)
  sum(width(hits))
}

#' Return "early" or "late" given strand and genomic position vs canonical
#' "early" = genomically upstream (lower coord on +, higher coord on -)
strand_shift_direction <- function(strand, position, canonical) {
  if (strand == "+") ifelse(position < canonical, "early", "late")
  else               ifelse(position > canonical, "early", "late")
}

#' Sample one shift distance from a truncated half-normal
#' Mean = shift_min (most shifts are small), SD = 50 bp, truncated to [shift_min, shift_max]
sample_shift <- function(shift_min, shift_max) {
  round(truncnorm::rtruncnorm(1, a = shift_min, b = shift_max,
                               mean = shift_min, sd = 50))
}

#' Format a novel transcript ID
make_novel_tx_id <- function(transcript_id, event_type) {
  paste0("NOVEL_", transcript_id, "_", event_type)
}

#' Consistent prefixed logging
log_info <- function(msg) cat(sprintf("[INFO] %s\n", msg))
log_warn <- function(msg) message(sprintf("[WARN] %s", msg))
