suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(here)
})
source(here::here("scripts/utils.R"))

#' Collapse a GTF to one transcript per gene
#'
#' Selects the transcript with the greatest total spliced length per gene.
#' Tiebreak: highest exon count. Excludes:
#'   - single-exon transcripts (< 3 exons)
#'   - transcripts lacking gene_id
#'   - chrM and MT
#'   - unassembled contigs (chromosome names containing "_")
#'
#' @param exons GRanges of exon features (type == "exon")
#' @return GRanges filtered to the winning transcript per gene
collapse_reference <- function(exons) {
  # Exclude chrM, MT, unassembled contigs, and entries without gene_id
  chr_names <- as.character(seqnames(exons))
  keep <- !chr_names %in% c("Y","chrM", "MT", "GL000213.1", "KI270711.1", "KI270721.1", "KI270727.1", "KI270734.1") &
          !grepl("_", chr_names) &
          !is.na(exons$gene_id)
  exons <- exons[keep]

  if (length(exons) == 0L) stop("No eligible exons after chromosome filtering.")

  # Compute per-transcript stats
  dt <- data.table(
    transcript_id = exons$transcript_id,
    gene_id       = exons$gene_id,
    width         = width(exons)
  )
  stats <- dt[, .(spliced_length = sum(width), exon_count = .N),
              by = .(transcript_id, gene_id)]

  # Require at least 3 exons (single-exon transcripts excluded per spec)
  stats <- stats[exon_count >= 3L]

  if (nrow(stats) == 0L) stop("No eligible transcripts with >= 3 exons.")

  # Pick best per gene: longest spliced length, then most exons as tiebreaker
  setorder(stats, gene_id, -spliced_length, -exon_count)
  best <- stats[, .SD[1L], by = gene_id]

  exons[exons$transcript_id %in% best$transcript_id]
}

# ── Snakemake glue ────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")

  log_info("Loading reference GTF...")
  exons <- load_gtf(snakemake@input[["gtf"]])
  log_info(sprintf("Loaded %d exons from %d transcripts",
                   length(exons), length(unique(exons$transcript_id))))

  collapsed <- collapse_reference(exons)
  log_info(sprintf("Collapsed to %d transcripts across %d genes",
                   length(unique(collapsed$transcript_id)),
                   length(unique(collapsed$gene_id))))

  # Re-import full GTF and filter to retained transcript + gene features
  full_gtf <- rtracklayer::import(snakemake@input[["gtf"]])
  keep_tx   <- full_gtf$transcript_id %in% collapsed$transcript_id
  keep_gene <- full_gtf$type == "gene" & full_gtf$gene_id %in% collapsed$gene_id
  full_gtf  <- full_gtf[keep_tx | keep_gene]

  rtracklayer::export(full_gtf, snakemake@output[[1]], format = "gtf")
  log_info("Done.")

  sink(); sink(type = "message")
  close(log_con)
}
