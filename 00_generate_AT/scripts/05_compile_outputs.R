suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(here)
  library(tidyverse)
})

source(here::here("scripts/utils.R"))

# ── Snakemake glue ────────────────────────────────────────────────────────────
# No testable pure-R logic here — this script is pure I/O assembly.
# Validated by the end-to-end integration test (Task 9).
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")

  collapsed_gtf <- rtracklayer::import(snakemake@input[["collapsed_gtf"]])
  novel_gtf     <- rtracklayer::import(snakemake@input[["novel_gtf"]])
  reference_gtf <- rtracklayer::import(snakemake@input[["reference_gtf"]])

  # novel_transcripts.gtf — novel isoforms only
  rtracklayer::export(novel_gtf, snakemake@output[["novel_transcripts"]], format = "gtf")
  log_info("Written novel_transcripts.gtf")

  # novel_transcript_ids.tsv — single-column, header "transcript_id"
  novel_ids <- unique(novel_gtf$transcript_id)
  write.table(data.frame(transcript_id = novel_ids),
              snakemake@output[["novel_ids"]],
              sep = "\t", quote = FALSE, row.names = FALSE)
  # strip everything but the original transcript ID of the finalnovel transcripts to
  # get original_transcript_ids.tsv — single-column, header "transcript_id"
  og_ids <- unique(novel_gtf$transcript_id) %>% 
    gsub("^NOVEL_", "", .) %>% 
    sub("_.*$", "", .)
  write.table(data.frame(transcript_id = og_ids),
              snakemake@output[["og_ids"]],
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_info("Written novel_transcript_ids.tsv and original_transcript_ids.tsv")

  og_gtf <- collapsed_gtf[collapsed_gtf$transcript_id %in% og_ids]
  
  # affected_gene_ids.tsv — single-column, header "gene_id"
  affected_genes <- unique(novel_gtf$gene_id)
  write.table(data.frame(gene_id = affected_genes),
              snakemake@output[["affected_genes"]],
              sep = "\t", quote = FALSE, row.names = FALSE)

  log_info(sprintf("Done. %d novel transcripts across %d genes.",
                   length(novel_ids), length(affected_genes)))

  # remove genes from reference which have og_ids to avoid duplicates in augmented GTF
  regular_gtf <- reference_gtf[!reference_gtf$gene_id %in% affected_genes]

  # augmented_transcriptome.gtf — novel + reference
  augmented <- c(regular_gtf, og_gtf, novel_gtf)
  rtracklayer::export(augmented, snakemake@output[["augmented"]], format = "gtf")
  log_info("Written augmented_transcriptome.gtf")
  
  control_gtf <- c(regular_gtf, og_gtf)
  rtracklayer::export(control_gtf, snakemake@output[["control"]], format = "gtf")
  log_info("Written control_transcriptome.gtf")

  sink(); sink(type = "message")
  close(log_con)
}
