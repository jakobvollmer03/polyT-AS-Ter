# function for bidirectional filtering between a strintie gtf and a set of differentially spliced junctions. 
# Used to compare accuracy current SpliCeAT pipeline against single DS detection tools on the transcript level 


# Libraries
suppressMessages({
    library(rtracklayer)
    library(GenomicFeatures)
    library(txdbmaker)
    library(data.table)
    library(dplyr)
})

tool <- "granges" # "majiq" # "leafcutter" #
consensus_events    <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/masterlists/consensus_granges_LSV.tsv"
stringtie_gtf       <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/stringtie_assemblies/stringtie_assembly.gtf"
reference_gtf       <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r0_get_ref/Mus_musculus/Mus_musculus_GRCm39_115_chr.gtf"
consensus_GTF_LSVs  <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/masterlists"
novel_gtf           <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/granges_stringtie_novel_transcripts.gtf"

# Remove ensembl version numbers and clean metadata columns for matching with LSV file
clean_ids <- function(x) {
    x <- gsub("transcript:|gene:", "", as.character(x))
    x <- ifelse(grepl("^MSTRG", x), x, sub("\\..*$", "", x))
    return(x) }

# 1. TRUTH SET (Ref + ds_detection consensus)
cat(sprintf("TRUTH SET: Processing Reference GTF and Granges Consensus File...\n"))
official_ref <- rtracklayer::import(reference_gtf)
#official_ref <- official_ref[official_ref$transcript_biotype == "protein_coding"]
ref_txdb <- txdbmaker::makeTxDbFromGRanges(official_ref)
ref_introns <- GenomicFeatures::intronsByTranscript(ref_txdb) %>% unlist() %>% unique()

consensus_matrix <- fread(
    consensus_events,
    header     = TRUE,
    fill       = TRUE,
    sep        = "\t",
    na.strings = character(0)   # nothing is auto-converted to NA; "NA" stays as the string "NA"
) %>% as.data.frame()
use_granges = TRUE
if (use_granges) {
  # Event Level Consensus
  cat(sprintf("TRUTH SET: Collapsing feature types from Granges Consensus File...\n"))
  collapse_feature_type <- function(m, w, l) {
    types <- c(as.character(m), as.character(w), as.character(l))
    types <- gsub("\\[|\\]", "", types)
    if (any(grepl("intron_retention|IR_sign_match", types, ignore.case = TRUE), na.rm = TRUE)) {
      return("intron_retention") }
    if (any(grepl("splice_junction|boundary match", types, ignore.case = TRUE), na.rm = TRUE)) {
      return("splice_junction")}
    return(NA_character_)
  }
  consensus_matrix <- consensus_matrix %>%
    rowwise() %>%
    mutate(feature_type = collapse_feature_type(
      majiq_feature_type,
      whippet_feature_type,
      leafcutter_feature_type
    )) %>%
    ungroup()
} else {
  cat(sprintf("TRUTH SET: Using gene-level consensus feature types directly...\n"))
  consensus_matrix <- consensus_matrix %>%
    rename(feature_type = feature_type)
}
consensus_matrix <- consensus_matrix %>%
  select(chr, strand, start, end, gene_id, gene_name, feature_type) %>%
  filter(!is.na(feature_type))

cat("TRUTH SET: Extracting ranges from Granges Consensus File...\n")
consensus_ranges <- GRanges(
    seqnames = consensus_matrix$chr,
    ranges = IRanges(start = consensus_matrix$start, end = consensus_matrix$end),
    strand = consensus_matrix$strand,
    feature_type = consensus_matrix$feature_type
)

sj_consensus <- consensus_ranges[consensus_ranges$feature_type == "splice_junction"]
cat("TRUTH SET: Splice Junctions Extracted.\n")
ir_consensus <- consensus_ranges[consensus_ranges$feature_type == "intron_retention"]
cat("TRUTH SET: Intron Retention events Extracted.\n")

# 2. TEST SET (StringTie)
cat("STRINGTIE: Building TxDb and extracting introns...\n")
gtf_raw <- rtracklayer::import(stringtie_gtf)
gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]

txdb_st <- txdbmaker::makeTxDbFromGRanges(gtf_clean)
st_introns_grl <- GenomicFeatures::intronsByTranscript(txdb_st, use.names=TRUE)
st_introns_flat <- unlist(st_introns_grl)

all_st_tx_names <- names(st_introns_grl)
novel_st_tx_names <- all_st_tx_names[grepl("^MSTRG", all_st_tx_names)]

st_exons_grl <- GenomicFeatures::exonsBy(txdb_st, by="tx", use.names=TRUE)
st_exons_flat <- unlist(st_exons_grl) 

# 3. FILTERING
cat("FILTER: Validating Splice Junctions and Intron Retentions...\n")

# 3.1. Remove transcripts with stringtie-only novel events
# 3.1.1 Event is splice junction
safe_SJ_ranges <- unique(c(ref_introns, sj_consensus)) #combined truth set
matches_safe <- findOverlaps(st_introns_flat, safe_SJ_ranges, type="equal", maxgap=1)

unsupported_sj_indices <- setdiff(seq_along(st_introns_flat), queryHits(matches_safe)) # stringtie SJ minus truth set SJ
unsupported_sj_tx <- unique(names(st_introns_flat[unsupported_sj_indices]))
cat(sprintf("FILTER: Splice Junction Validation Complete. Removed %d transcripts with unsupported splice junctions.\n", length(unsupported_sj_tx)))

# §3.1.2 — IR detection and validation
# Primary signal: does any StringTie exon fully contain a MAJIQ IR event?
hits_ir  <- findOverlaps(st_exons_flat, ir_consensus, type = "any", ignore.strand = FALSE)
exon_ir  <- st_exons_flat[queryHits(hits_ir)]
ir_hit   <- ir_consensus[subjectHits(hits_ir)]

exon_contains_ir   <- start(exon_ir) <= start(ir_hit) + 5 &
                      end(exon_ir)   >= end(ir_hit)   - 5
tx_with_ir_support <- unique(names(exon_ir)[exon_contains_ir])

# Unsupported IR: exon spans a reference intron but has NO MAJIQ backing
hits_ref          <- findOverlaps(st_exons_flat, ref_introns, type = "any", ignore.strand = FALSE)
exon_ref          <- st_exons_flat[queryHits(hits_ref)]
intron_ref        <- ref_introns[subjectHits(hits_ref)]

exon_spans_intron       <- start(exon_ref) <= start(intron_ref) + 5 &
                           end(exon_ref)   >= end(intron_ref)   - 5
tx_spanning_ref_introns <- unique(names(exon_ref)[exon_spans_intron])

unsupported_ir_tx <- setdiff(tx_spanning_ref_introns, tx_with_ir_support)
cat(sprintf("FILTER: Intron Retention Validation Complete. Removed %d transcripts with unsupported intron retention events.\n", length(unsupported_ir_tx)))

# §3.2.1 — novel SJ check (unchanged)
matches_sj_lsv     <- findOverlaps(st_introns_flat, sj_consensus, type="equal", maxgap=1)
tx_with_sj_support <- unique(names(st_introns_flat[queryHits(matches_sj_lsv)]))

# §3.2.2 — novel IR check
hits_ir_lsv        <- findOverlaps(st_exons_flat, ir_consensus, type="any", ignore.strand=FALSE)
exon_ir_lsv        <- st_exons_flat[queryHits(hits_ir_lsv)]
ir_lsv_hit         <- ir_consensus[subjectHits(hits_ir_lsv)]
contains_ir_lsv    <- start(exon_ir_lsv) <= start(ir_lsv_hit) + 5 &
                      end(exon_ir_lsv)   >= end(ir_lsv_hit)   - 5
tx_with_ir_support <- unique(names(exon_ir_lsv)[contains_ir_lsv])

# Combine — this line was missing
supported_by_truth <- unique(c(tx_with_sj_support, tx_with_ir_support))

unsupported_tx_names <- unique(c(unsupported_sj_tx, unsupported_ir_tx))
cat(sprintf("FILTER: %d transcripts contain unvalidated events. They have been removed.\n", length(unsupported_tx_names)))
# 3.3 Combined filters to identify valid transcripts
valid_tx_names <- intersect(supported_by_truth, setdiff(novel_st_tx_names, unsupported_tx_names))
cat(sprintf("FILTER: Overall Filtering Complete. %d transcripts contains at least one validated novel event.\n", length(valid_tx_names)))

# 4. EXPORT
# 4.1. EXPORT CONSENSUS MATRIX
cat("EXPORT: Generating validated consensus matrix...\n")

#Subset features to the final "winning" novel transcripts
st_final_introns <- st_introns_flat[names(st_introns_flat) %in% valid_tx_names]
st_final_exons      <- st_exons_grl[names(st_exons_grl) %in% valid_tx_names]
st_final_exons_flat <- unlist(st_final_exons)

hits_ir_export <- findOverlaps(st_final_exons_flat, ir_consensus, type = "any", ignore.strand = FALSE)
exon_exp       <- st_final_exons_flat[queryHits(hits_ir_export)]
ir_exp         <- ir_consensus[subjectHits(hits_ir_export)]

contained      <- start(exon_exp) <= start(ir_exp) + 5 &
                  end(exon_exp)   >= end(ir_exp)   - 5
valid_sj_idx <- unique(queryHits(findOverlaps(sj_consensus, st_final_introns, type="equal", maxgap=1)))
valid_ir_idx   <- unique(subjectHits(hits_ir_export)[contained])  # index into ir_consensus → consensus_matrix

final_consensus_matrix <- rbind(
    consensus_matrix %>% filter(feature_type == "splice_junction") %>% slice(valid_sj_idx),
    consensus_matrix %>% filter(feature_type == "intron_retention") %>% slice(valid_ir_idx)
)

cat(sprintf("EXPORT: %d / %d %s events confirmed by novel StringTie transcripts.\n",
            nrow(final_consensus_matrix), nrow(consensus_matrix), tool))

consensus_outfile <- file.path(consensus_GTF_LSVs, paste0(tool, "_confirmed_events_full_ref.tsv"))
fwrite(final_consensus_matrix, consensus_outfile, sep="\t", na="NA", quote=FALSE)
cat(sprintf("EXPORT: Consensus Matrix Complete. File Saved at %s\n", consensus_GTF_LSVs))

# 4.2. EXPORT GTF
# 4.2.1. Novel transcripts only (filtered StringTie)
cat(sprintf("EXPORT: Generating validated consensus matrix (supported by StringTie)..."))
gtf_novel_only <- gtf_clean[mcols(gtf_clean)$transcript_id %in% valid_tx_names]
rtracklayer::export(gtf_novel_only, novel_gtf, format="gtf")
cat(sprintf("EXPORT: Novel-only GTF Construction Complete. File Saved at %s",novel_gtf))