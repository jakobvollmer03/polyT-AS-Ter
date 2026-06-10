suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(here)
})

source(here::here("scripts/utils.R"))

# ── Core generator ────────────────────────────────────────────────────────────

#' Generate one novel transcript from a candidate row.
#'
#' @param candidate  One-row data.frame from 02_select_candidates output
#' @param exons      GRanges of exon features (collapsed reference, exon rows only)
#' @return list with:
#'   $novel_tx_id  character
#'   $novel_exons  data.frame(chr, strand, gene_id, transcript_id, exon_start, exon_end, exon_number)
#'   $fasta_mod    data.frame(chr, pos1, pos2, base1, base2, transcript_id, site) or NULL (IR)
generate_one <- function(candidate, exons) {
  tx_id    <- candidate$transcript_id
  gene_id  <- candidate$gene_id
  etype    <- candidate$event_type
  strand   <- candidate$strand
  novel_id <- make_novel_tx_id(tx_id, etype)

  # FIX Bug A (partial): filter to exon-type features before subsetting.
  # If load_gtf returns mixed feature types, transcript/gene rows contaminate
  # tx_gr and produce spurious large-span "exons" that corrupt coordinate logic
  # and can produce NA intron lengths during candidate building (if not already
  # filtered in step 2). Belt-and-suspenders filter applied here too.
  exons_only <- if (!is.null(exons$type)) exons[exons$type == "exon"] else exons
  tx_gr <- exons_only[exons_only$transcript_id == tx_id]
  if (length(tx_gr) == 0L)
    stop(sprintf("generate_one: no exon features found for transcript '%s'. ",
                 "Check that the collapsed GTF contains exon-level rows and that ",
                 "transcript IDs match between the GTF and candidates table.", tx_id))

  tx_df    <- as.data.frame(tx_gr)
  tx_df    <- tx_df[order(tx_df$start), ]

  tx_exons <- data.frame(
    chr           = as.character(tx_df$seqnames),
    strand        = as.character(tx_df$strand),
    gene_id       = as.character(tx_df$gene_id),
    transcript_id = novel_id,
    exon_start    = as.integer(tx_df$start),
    exon_end      = as.integer(tx_df$end),
    stringsAsFactors = FALSE
  )

  # Sanity check: no NA coordinates allowed before we start transforming
  if (any(is.na(tx_exons$exon_start)) || any(is.na(tx_exons$exon_end)))
    stop(sprintf(
      "generate_one: NA coordinates in source exons for transcript '%s'. ",
      "This indicates a corrupt collapsed GTF or an exon-type filter failure.",
      tx_id))

  base_type <- sub("_ES_[12]$", "", etype)
  is_combined <- grepl("_ES_[12]$", etype)
  n_skip_combined <- if (grepl("_ES_2$", etype)) 2L else if (grepl("_ES_1$", etype)) 1L else 0L

  fasta_mod <- NULL

  # ── Apply alt boundary component ───────────────────────────────────────────
  if (base_type %in% c("alternative_donor_late",    "alternative_donor_early",
                        "alternative_acceptor_late", "alternative_acceptor_early")) {

    # FIX Bug A (root): validate shift values before use.
    # shift_max = NA occurs when the candidate table contains a row built from
    # a non-exon feature row (producing a nonsensical intron_len), which then
    # propagates as NA through sample_shift into coordinate arithmetic.
    if (is.na(candidate$shift_min) || is.na(candidate$shift_max))
      stop(sprintf(
        "generate_one: shift_min or shift_max is NA for candidate gene '%s' ",
        "event '%s'. This indicates a contaminated candidate row — likely a ",
        "non-exon feature row was not filtered out in step 2. Check that ",
        "load_gtf and 02_select_candidates.R both filter to type == 'exon'.",
        gene_id, etype))

    shift  <- sample_shift(candidate$shift_min, candidate$shift_max)

    # Guard: sample_shift must not return NA (e.g. from truncnorm on sd=0 edge case)
    if (is.na(shift) || is.nan(shift)) {
      message(sprintf(
        "  [WARN] sample_shift returned NA/NaN for shift_min=%d shift_max=%d ",
        "(gene %s, event %s). Falling back to shift_min.",
        as.integer(candidate$shift_min), as.integer(candidate$shift_max), gene_id, etype))
      shift <- as.integer(candidate$shift_min)
    }

    result   <- apply_alt_boundary(tx_exons, candidate, base_type, strand, shift, novel_id)
    tx_exons  <- result$exons
    fasta_mod <- result$fasta_mod
  }

  # ── Apply exon-skip component (covers pure ES and combined events) ──────────
  if (etype %in% c("1_ES", "2_ES") || is_combined) {
    n_skip <- if (etype == "1_ES") 1L else if (etype == "2_ES") 2L else n_skip_combined

    # For combined events the skip region is not stored in the candidate row
    # (the candidate was built for the alt boundary, not the skip).
    # We therefore select the skip target dynamically: for combined events,
    # pick a consecutive internal exon run that is NOT the alt-boundary exon.
    if (is_combined) {
      result <- apply_exon_skip_dynamic(tx_exons, candidate, n_skip, novel_id, gene_id)
    } else {
      result <- apply_exon_skip(tx_exons, candidate, n_skip, novel_id)
    }
    tx_exons <- result$exons

    # Accumulate fasta_mod: combined events produce two junction modification sets
    if (!is.null(result$fasta_mod)) {
      fasta_mod <- if (is.null(fasta_mod)) result$fasta_mod else rbind(fasta_mod, result$fasta_mod)
    }
  }

  if (etype == "intron_retention") {
    result   <- apply_ir(tx_exons, candidate, novel_id)
    tx_exons <- result$exons
    # IR creates no new splice junction — fasta_mod stays NULL
  }

  tx_exons$exon_number <- seq_len(nrow(tx_exons))

  # Final validation before handing off to exons_to_granges
  validate_novel_exons(tx_exons, novel_id)

  list(novel_tx_id = novel_id, novel_exons = tx_exons, fasta_mod = fasta_mod)
}

# ── Event transformers ────────────────────────────────────────────────────────

#' Apply an alternative donor or acceptor boundary shift.
#'
#' "late"  = shift into the intron (exon gets longer).
#' "early" = shift into the exon  (exon gets shorter).
#' Directions are strand-aware:
#'   + strand donor   : exon_end   moves right (late) / left  (early)
#'   + strand acceptor: exon_start moves left  (late) / right (early)
#'   - strand donor   : exon_start moves left  (late) / right (early)
#'   - strand acceptor: exon_end   moves right (late) / left  (early)
apply_alt_boundary <- function(tx_exons, candidate, base_type, strand, shift, novel_id) {
  te_start <- as.integer(candidate$target_exon_start)
  te_end   <- as.integer(candidate$target_exon_end)

  # FIX Bug D: validate target exon lookup — integer(0) causes silent no-op
  i <- which(tx_exons$exon_start == te_start & tx_exons$exon_end == te_end)
  if (length(i) == 0L)
    stop(sprintf(
      "apply_alt_boundary: target exon [%d, %d] not found in transcript '%s'. ",
      "This can occur if exon coordinates were altered by a 0-based/1-based ",
      "conversion during GTF export/import, or if non-exon rows are present. ",
      "Check coordinate system consistency throughout the pipeline.",
      te_start, te_end, candidate$transcript_id))
  if (length(i) > 1L)
    stop(sprintf(
      "apply_alt_boundary: target exon [%d, %d] matches %d rows in transcript '%s'. ",
      "Exon coordinates must be unique within a transcript.",
      te_start, te_end, length(i), candidate$transcript_id))

  late  <- grepl("late",  base_type)
  donor <- grepl("donor", base_type)

  new_start <- te_start
  new_end   <- te_end

  if (strand == "+") {
    if ( donor &&  late) new_end   <- te_end   + shift   # right into intron
    if ( donor && !late) new_end   <- te_end   - shift   # left  into exon
    if (!donor &&  late) new_start <- te_start - shift   # left  into intron
    if (!donor && !late) new_start <- te_start + shift   # right into exon
  } else {
    # Minus strand: donor is at low coord (exon_start), acceptor at high (exon_end)
    if ( donor &&  late) new_start <- te_start - shift   # left  into intron
    if ( donor && !late) new_start <- te_start + shift   # right into exon
    if (!donor &&  late) new_end   <- te_end   + shift   # right into intron
    if (!donor && !late) new_end   <- te_end   - shift   # left  into exon
  }

  # Derive novel intron boundaries for FASTA modification
  if (strand == "+") {
    if (donor) {
      I_start <- new_end + 1L
      I_end   <- as.integer(candidate$target_intron_end)
    } else {
      I_start <- as.integer(candidate$target_intron_start)
      I_end   <- new_start - 1L
    }
  } else {
    if (donor) {
      I_start <- as.integer(candidate$target_intron_start)
      I_end   <- new_start - 1L
    } else {
      I_start <- new_end + 1L
      I_end   <- as.integer(candidate$target_intron_end)
    }
  }

  if (I_end < I_start)
    stop(sprintf(
      "apply_alt_boundary: derived intron [%d, %d] is invalid (end < start) for ",
      "gene '%s', event '%s', shift=%d. Check coordinate logic for strand '%s'.",
      I_start, I_end, candidate$gene_id, candidate$event_type, shift, strand))

  tx_exons[i, "exon_start"] <- new_start
  tx_exons[i, "exon_end"]   <- new_end

  fmod <- fasta_mod_for_junction(strand, I_start, I_end, novel_id, candidate$chr)
  list(exons = tx_exons, fasta_mod = fmod)
}

#' Remove n_skip consecutive internal exons (exon skipping event) using
#' the skip region stored in the candidate row (pure ES events only).
apply_exon_skip <- function(tx_exons, candidate, n_skip, novel_id) {
  skip_start <- as.integer(candidate$target_exon_start)  # upstream_exon_end + 1
  upstream_hits <- which(tx_exons$exon_end < skip_start)
  if (length(upstream_hits) == 0L)
    stop(sprintf("apply_exon_skip: no upstream exon found for skip_start=%d", skip_start))
  i_upstream  <- max(upstream_hits)
  remove_idx  <- (i_upstream + 1L):(i_upstream + n_skip)
  remove_idx  <- remove_idx[remove_idx <= nrow(tx_exons)]
  if (length(remove_idx) == 0L)
    stop(sprintf(
      "apply_exon_skip: no exons to remove (i_upstream=%d, n_skip=%d, n_exons=%d)",
      i_upstream, n_skip, nrow(tx_exons)))
  tx_exons <- tx_exons[-remove_idx, ]
  if (i_upstream + 1L > nrow(tx_exons))
    stop(sprintf(
      "apply_exon_skip: no downstream exon after removal (i_upstream=%d, remaining=%d)",
      i_upstream, nrow(tx_exons)))
  I_start <- tx_exons[i_upstream,      "exon_end"]   + 1L
  I_end   <- tx_exons[i_upstream + 1L, "exon_start"] - 1L
  fmod    <- fasta_mod_for_junction(tx_exons$strand[1L], I_start, I_end,
                                    novel_id, tx_exons$chr[1L])
  list(exons = tx_exons, fasta_mod = fmod)
}

#' For combined events: pick an internal exon run to skip that does NOT
#' include the exon that was already modified by apply_alt_boundary.
#' The modified exon is identified by comparing current coords to the
#' candidate's target coords (which no longer match after the boundary shift).
apply_exon_skip_dynamic <- function(tx_exons, candidate, n_skip, novel_id, gene_id) {
  # The alt-boundary exon has already been modified: its coords differ from
  # candidate$target_exon_start / target_exon_end.
  # Identify which row was modified: it will NOT match the original target coords.
  te_start <- as.integer(candidate$target_exon_start)
  te_end   <- as.integer(candidate$target_exon_end)
  modified_row <- which(!(tx_exons$exon_start == te_start & tx_exons$exon_end == te_end) &
                        seq_len(nrow(tx_exons)) %in% c(
                          which(tx_exons$exon_start == te_start),
                          which(tx_exons$exon_end   == te_end)))
  # Simpler: find internal exon runs of length n_skip that do not include the
  # modified exon. We identify the modified exon as the one closest to the
  # original target_exon_start (its start or end will differ by exactly `shift`).
  n <- nrow(tx_exons)
  # Collect candidate skip windows: runs of n_skip consecutive internal exons
  # (index 2 to n-1) that exclude whatever row shares the original target position
  skip_windows <- list()
  for (i_up in 1L:(n - n_skip - 1L)) {
    run <- (i_up + 1L):(i_up + n_skip)
    # Exclude windows that contain the alt-boundary exon
    # Alt-boundary exon is identified as the row whose ORIGINAL coords match
    # target_exon_start/end — i.e., the row nearest those coords.
    # Since that row is now modified, check whether any run row's CURRENT
    # coords are close to the original target (within 300 bp as a proxy).
    contains_modified <- any(
      abs(tx_exons$exon_start[run] - te_start) < 300L |
      abs(tx_exons$exon_end[run]   - te_end)   < 300L
    )
    if (!contains_modified) skip_windows[[length(skip_windows) + 1L]] <- i_up
  }
  if (length(skip_windows) == 0L)
    stop(sprintf(
      "apply_exon_skip_dynamic: no valid exon skip window found for ",
      "combined event in gene '%s'. Transcript may not have enough exons.",
      gene_id))
  # Pick a random valid window
  i_upstream  <- unlist(sample(skip_windows, 1L))
  remove_idx  <- (i_upstream + 1L):(i_upstream + n_skip)
  tx_exons    <- tx_exons[-remove_idx, ]
  I_start     <- tx_exons[i_upstream,      "exon_end"]   + 1L
  I_end       <- tx_exons[i_upstream + 1L, "exon_start"] - 1L
  fmod        <- fasta_mod_for_junction(tx_exons$strand[1L], I_start, I_end,
                                        novel_id, tx_exons$chr[1L])
  list(exons = tx_exons, fasta_mod = fmod)
}

#' Merge two adjacent exons by retaining the intron between them.
apply_ir <- function(tx_exons, candidate, novel_id) {
  I_start  <- as.integer(candidate$target_intron_start)
  hits     <- which(tx_exons$exon_end < I_start)
  if (length(hits) == 0L)
    stop(sprintf("apply_ir: no left exon found for intron_start=%d", I_start))
  i_left   <- max(hits)
  i_right  <- i_left + 1L
  if (i_right > nrow(tx_exons))
    stop(sprintf("apply_ir: no right exon found (i_left=%d, n_exons=%d)",
                 i_left, nrow(tx_exons)))
  tx_exons[i_left, "exon_end"] <- tx_exons[i_right, "exon_end"]
  tx_exons <- tx_exons[-i_right, ]
  list(exons = tx_exons, fasta_mod = NULL)
}

# ── FASTA modification helper ─────────────────────────────────────────────────

#' Compute the genomic positions that need GT-AG substitution at a novel junction.
#'
#' Returns a 2-row data.frame (donor row + acceptor row).
#' Positions are 1-based genomic coordinates.
#'
#' GT-AG rules (from CLAUDE.md):
#'   + strand donor   : genome[I_start]   = G, genome[I_start+1] = T
#'   + strand acceptor: genome[I_end-1]   = A, genome[I_end]     = G
#'   - strand donor   : genome[I_end-1]   = A, genome[I_end]     = C  (CT on +, GT on -)
#'   - strand acceptor: genome[I_start]   = C, genome[I_start+1] = T  (CT on +, AG on -)
fasta_mod_for_junction <- function(strand, I_start, I_end, tx_id, chr) {
  if (strand == "+") {
    donor_mod    <- data.frame(chr = chr, pos1 = I_start,     pos2 = I_start + 1L,
                               base1 = "G", base2 = "T",
                               transcript_id = tx_id, site = "donor",
                               stringsAsFactors = FALSE)
    acceptor_mod <- data.frame(chr = chr, pos1 = I_end - 1L,  pos2 = I_end,
                               base1 = "A", base2 = "G",
                               transcript_id = tx_id, site = "acceptor",
                               stringsAsFactors = FALSE)
  } else {
    donor_mod    <- data.frame(chr = chr, pos1 = I_end - 1L,  pos2 = I_end,
                               base1 = "A", base2 = "C",
                               transcript_id = tx_id, site = "donor",
                               stringsAsFactors = FALSE)
    acceptor_mod <- data.frame(chr = chr, pos1 = I_start,     pos2 = I_start + 1L,
                               base1 = "C", base2 = "T",
                               transcript_id = tx_id, site = "acceptor",
                               stringsAsFactors = FALSE)
  }
  rbind(donor_mod, acceptor_mod)
}

# ── Validation helper ─────────────────────────────────────────────────────────

#' Abort with a clear message if any exon coordinate is NA or invalid.
validate_novel_exons <- function(exons_df, tx_id) {
  if (any(is.na(exons_df$exon_start)) || any(is.na(exons_df$exon_end)))
    stop(sprintf("validate_novel_exons: NA coordinates in novel exons for '%s'", tx_id))
  if (any(exons_df$exon_end < exons_df$exon_start))
    stop(sprintf(
      "validate_novel_exons: exon_end < exon_start in novel exons for '%s'. ",
      "Rows: %s",
      tx_id,
      paste(which(exons_df$exon_end < exons_df$exon_start), collapse = ",")))
}

# ── GTF export helper ─────────────────────────────────────────────────────────

#' Convert a data.frame of novel exons to GRanges suitable for GTF export.
exons_to_granges <- function(exons_df, gene_id, tx_id) {
  gr <- GRanges(
    seqnames = exons_df$chr,
    ranges   = IRanges(start = exons_df$exon_start, end = exons_df$exon_end),
    strand   = exons_df$strand
  )
  gr$type          <- "exon"
  gr$gene_id       <- gene_id
  gr$transcript_id <- tx_id
  gr$exon_number   <- exons_df$exon_number
  gr
}

# ── Snakemake glue ────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")

  set.seed(snakemake@params[["seed"]])

  exons <- load_gtf(snakemake@input[["collapsed_gtf"]])

  # Diagnostic: report feature type breakdown so missing exon filter is visible
  message(sprintf("[INFO] Loaded collapsed GTF: %d features", length(exons)))
  if (!is.null(exons$type)) {
    type_tab <- sort(table(as.character(exons$type)), decreasing = TRUE)
    message("[INFO] Feature types: ",
            paste(names(type_tab), type_tab, sep = "=", collapse = ", "))
    n_before <- length(exons)
    exons    <- exons[exons$type == "exon"]
    message(sprintf("[INFO] After exon-type filter: %d features retained", length(exons)))
    if (length(exons) == 0L)
      stop("[FATAL] No exon-type features remain after filtering. ",
           "Check that the collapsed GTF contains rows with type == 'exon'.")
  } else {
    message("[WARN] No 'type' column in collapsed GTF — assuming all rows are exons")
  }

  candidates <- read.table(snakemake@input[["candidates"]], header = TRUE,
                           sep = "\t", stringsAsFactors = FALSE)
  message(sprintf("[INFO] Processing %d candidates", nrow(candidates)))

  novel_gr_list <- list()
  fasta_mods    <- list()

  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, ]
    log_info(sprintf("Generating %s for gene %s", cand$event_type, cand$gene_id))
    result <- tryCatch(
      generate_one(cand, exons),
      error = function(e) {
        log_warn(sprintf("FAILED %s gene %s: %s",
                         cand$event_type, cand$gene_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(result)) next

    novel_gr_list[[length(novel_gr_list) + 1L]] <-
      exons_to_granges(result$novel_exons, cand$gene_id, result$novel_tx_id)

    # FIX Bug C: only append non-NULL fasta_mod entries — IR events return NULL
    if (!is.null(result$fasta_mod))
      fasta_mods[[length(fasta_mods) + 1L]] <- result$fasta_mod
  }

  if (length(novel_gr_list) == 0L)
    stop("[FATAL] No novel transcripts were generated successfully. ",
         "Review FAILED lines above.")

  # FIX Bug C: do.call(rbind, list()) on a list that contained NULLs was
  # previously able to error or silently drop rows depending on R version.
  # NULLs are now excluded above; guard remains for safety.
  novel_gr <- do.call(c, novel_gr_list)
  rtracklayer::export(novel_gr, snakemake@output[["novel_gtf"]], format = "gtf")

  if (length(fasta_mods) > 0L) {
    mods_df <- do.call(rbind, fasta_mods)
    write.table(mods_df, snakemake@output[["fasta_mods"]],
                sep = "\t", quote = FALSE, row.names = FALSE)
    log_info(sprintf("Written %d novel transcripts, %d FASTA modification rows",
                     length(novel_gr_list), nrow(mods_df)))
  } else {
    # All events were IR — write header-only TSV so downstream step doesn't fail
    mods_df <- data.frame(chr = character(), pos1 = integer(), pos2 = integer(),
                          base1 = character(), base2 = character(),
                          transcript_id = character(), site = character())
    write.table(mods_df, snakemake@output[["fasta_mods"]],
                sep = "\t", quote = FALSE, row.names = FALSE)
    log_info(sprintf("Written %d novel transcripts, 0 FASTA modification rows (IR-only run)",
                     length(novel_gr_list)))
  }

  sink(); sink(type = "message")
  close(log_con)
}