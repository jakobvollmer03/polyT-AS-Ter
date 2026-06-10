suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(here)
})
source(here::here("scripts/utils.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# Minimum exon counts per event type (base types only)
EVENT_MIN_EXONS <- list(
  alternative_donor_late      = 3L,
  alternative_donor_early     = 3L,
  alternative_acceptor_late   = 3L,
  alternative_acceptor_early  = 3L,
  intron_retention            = 3L,
  "1_ES"                      = 4L,   # N_skipped + 3
  "2_ES"                      = 5L
)

# For combined types (e.g. alternative_donor_late_ES_1), take max of components
combined_min_exons <- function(event_type) {
  if (grepl("_ES_2$", event_type)) return(5L)
  if (grepl("_ES_1$", event_type)) return(4L)
  base <- sub("_ES_[12]$", "", event_type)
  EVENT_MIN_EXONS[[base]] %||% 3L
}

# ── Builder helpers ───────────────────────────────────────────────────────────

build_alt_boundary_candidates <- function(dt, base_type, event_type) {
  late  <- grepl("late",  base_type)
  donor <- grepl("donor", base_type)

  results <- list()
  for (tx_id in unique(dt$transcript_id)) {
    tx <- dt[transcript_id == tx_id][order(exon_start)]
    n  <- nrow(tx)
    if (n < 3L) next  # FIX Bug 4: guard against 2:(n-1) iterating backwards

    # Iterate only over internal exons (not first or last in genomic order)
    for (i in 2L:(n - 1L)) {
      exon <- tx[i]

      if (donor) {
        # Donor = the exon end abutting the downstream intron.
        # On + strand the intron is to the RIGHT (higher coords).
        # On - strand the intron is to the LEFT  (lower coords).
        if (exon$strand == "+") {
          intron_start <- exon$exon_end + 1L
          intron_end   <- tx[i + 1L]$exon_start - 1L
        } else {
          intron_start <- tx[i - 1L]$exon_end + 1L
          intron_end   <- exon$exon_start - 1L
        }
      } else {
        # Acceptor = the exon end abutting the upstream intron.
        # On + strand the intron is to the LEFT  (lower coords).
        # On - strand the intron is to the RIGHT (higher coords).
        if (exon$strand == "+") {
          intron_start <- tx[i - 1L]$exon_end + 1L
          intron_end   <- exon$exon_start - 1L
        } else {
          intron_start <- exon$exon_end + 1L
          intron_end   <- tx[i + 1L]$exon_start - 1L
        }
      }

      intron_len <- intron_end - intron_start + 1L
      exon_len   <- exon$exon_end - exon$exon_start + 1L

      # Sanity check: coordinates must be positive
      if (intron_len <= 0L || exon_len <= 0L) next

      if (late) {
        # Shift into intron → exon gets longer.
        # Intron must be >= 150 bp: 50 bp minimum shift + 100 bp remaining intron.
        if (intron_len < 150L) next
        s_min <- 50L
        s_max <- as.integer(min(250L, intron_len - 100L))
      } else {
        # Shift into exon → exon gets shorter.
        # Exon must be >= 100 bp: 50 bp minimum shift + 50 bp remaining exon.
        if (exon_len < 100L) next
        s_min <- 50L
        s_max <- as.integer(min(250L, exon_len - 50L))
      }

      if (s_max < s_min) next  # No valid shift range — skip

      results[[length(results) + 1L]] <- data.frame(
        gene_id             = exon$gene_id,
        transcript_id       = tx_id,
        event_type          = event_type,
        chr                 = exon$chr,
        strand              = exon$strand,
        target_exon_start   = exon$exon_start,
        target_exon_end     = exon$exon_end,
        target_intron_start = intron_start,
        target_intron_end   = intron_end,
        shift_min           = s_min,
        shift_max           = s_max,
        stringsAsFactors    = FALSE
      )
    }
  }
  if (length(results) == 0L) return(data.frame())
  do.call(rbind, results)
}

build_ir_candidates <- function(dt) {
  results <- list()
  for (tx_id in unique(dt$transcript_id)) {
    tx <- dt[transcript_id == tx_id][order(exon_start)]
    n  <- nrow(tx)
    if (n < 3L) next
    for (i in 1L:(n - 1L)) {
      intron_start <- tx[i]$exon_end + 1L
      intron_end   <- tx[i + 1L]$exon_start - 1L
      intron_len   <- intron_end - intron_start + 1L
      if (intron_len < 200L) next
      results[[length(results) + 1L]] <- data.frame(
        gene_id             = tx[i]$gene_id,
        transcript_id       = tx_id,
        event_type          = "intron_retention",
        chr                 = tx[i]$chr,
        strand              = tx[i]$strand,
        target_exon_start   = tx[i]$exon_start,
        target_exon_end     = tx[i + 1L]$exon_end,
        target_intron_start = intron_start,
        target_intron_end   = intron_end,
        shift_min           = NA_integer_,
        shift_max           = NA_integer_,
        stringsAsFactors    = FALSE
      )
    }
  }
  if (length(results) == 0L) return(data.frame())
  do.call(rbind, results)
}

build_es_candidates <- function(dt, event_type) {
  n_skip <- if (event_type == "1_ES") 1L else 2L
  results <- list()
  for (tx_id in unique(dt$transcript_id)) {
    tx <- dt[transcript_id == tx_id][order(exon_start)]
    n  <- nrow(tx)
    if (n < n_skip + 3L) next
    for (i in 1L:(n - n_skip - 1L)) {
      upstream_end     <- tx[i]$exon_end
      downstream_start <- tx[i + n_skip + 1L]$exon_start
      results[[length(results) + 1L]] <- data.frame(
        gene_id             = tx[i]$gene_id,
        transcript_id       = tx_id,
        event_type          = event_type,
        chr                 = tx[i]$chr,
        strand              = tx[i]$strand,
        target_exon_start   = upstream_end + 1L,
        target_exon_end     = downstream_start - 1L,
        target_intron_start = upstream_end + 1L,
        target_intron_end   = downstream_start - 1L,
        shift_min           = NA_integer_,
        shift_max           = NA_integer_,
        stringsAsFactors    = FALSE
      )
    }
  }
  if (length(results) == 0L) return(data.frame())
  do.call(rbind, results)
}

# ── Overlap filtering ─────────────────────────────────────────────────────────

#' Identify transcripts with no overlapping exons from other transcripts.
#'
#' Checks each transcript to see if its exons overlap (on the same chromosome)
#' with exons from any other transcript. Returns only transcript IDs that have
#' no such overlaps.
#'
#' @param exons  GRanges of exon features from the collapsed reference
#' @return Character vector of transcript IDs with no overlapping exons
get_non_overlapping_transcripts <- function(exons) {
  # Filter to exon features only to ensure consistent data
  if (!is.null(exons$type)) {
    exons <- exons[exons$type == "exon"]
  }
  if (length(exons) == 0L) return(character(0))

  # Build data.table for efficient overlap detection
  exons_dt <- data.table(
    chr           = as.character(seqnames(exons)),
    start         = start(exons),
    end           = end(exons),
    transcript_id = as.character(exons$transcript_id)
  )

  # Track transcripts that have overlaps with other transcripts
  overlapping_tx <- character(0)

  for (tx_id in unique(exons_dt$transcript_id)) {
    tx_exons   <- exons_dt[transcript_id == tx_id]
    other_exons <- exons_dt[transcript_id != tx_id]

    # Check if any exon from this transcript overlaps with other transcripts
    has_overlap <- FALSE
    for (i in seq_len(nrow(tx_exons))) {
      hit <- other_exons[
        chr == tx_exons[i, chr] &
        start <= tx_exons[i, end] &
        end >= tx_exons[i, start]
      ]
      if (nrow(hit) > 0L) {
        has_overlap <- TRUE
        break
      }
    }

    if (has_overlap) {
      overlapping_tx <- c(overlapping_tx, tx_id)
    }
  }

  # Return transcript IDs without overlaps
  all_tx <- unique(exons_dt$transcript_id)
  all_tx[!(all_tx %in% overlapping_tx)]
}

# ── Public API ────────────────────────────────────────────────────────────────

#' Return eligible candidate rows for a given event type.
#'
#' Each row represents one (transcript, target site) pair with pre-computed
#' shift_min/shift_max bounds.  The 500 bp spliced-length filter and the
#' per-type minimum exon count are applied here.
#'
#' @param exons  GRanges of features from the collapsed reference.
#'               Must include only exon-type rows — this function will
#'               attempt to filter automatically if a 'type' column is present.
#' @param event_type  Character scalar naming the event type
#' @return data.frame (zero-row if nothing eligible)
eligible_for_type <- function(exons, event_type) {

  # FIX Bug 2: Ensure we work only with exon-level features.
  # If the GRanges contains mixed feature types (gene/transcript/exon/CDS/UTR),
  # non-exon rows corrupt exon counts and spliced-length calculations.
  if (!is.null(exons$type)) {
    n_before <- length(exons)
    exons    <- exons[exons$type == "exon"]
    n_after  <- length(exons)
    if (n_after == 0L) {
      message(sprintf(
        "  [DIAG] eligible_for_type('%s'): filtered from %d to 0 features — ",
        event_type, n_before),
        "the GTF contains no features with type == 'exon'. ",
        "Check that 01_collapse_reference.R exports both transcript and exon ",
        "level rows, or that load_gtf() does not pre-filter incorrectly."
      )
      return(data.frame())
    }
    message(sprintf(
      "  [DIAG] eligible_for_type('%s'): retained %d / %d exon features",
      event_type, n_after, n_before
    ))
  }

  # FIX Bug 3: removed dead exon_number column — it is unused downstream and
  # as.integer(NULL) returns integer(0) which causes a data.table length error
  # if the collapsed GTF lacks this attribute.
  dt <- data.table(
    chr           = as.character(seqnames(exons)),
    strand        = as.character(strand(exons)),
    gene_id       = as.character(exons$gene_id),
    transcript_id = as.character(exons$transcript_id),
    exon_start    = start(exons),
    exon_end      = end(exons)
  )

  # Drop rows with missing identifiers — these indicate non-exon features that
  # slipped through without a type column (e.g. gene rows with NA transcript_id)
  dt <- dt[!is.na(transcript_id) & !is.na(gene_id) & transcript_id != "" & gene_id != ""]

  message(sprintf(
    "  [DIAG] eligible_for_type('%s'): %d exon rows, %d unique transcripts",
    event_type, nrow(dt), length(unique(dt$transcript_id))
  ))

  if (nrow(dt) == 0L) return(data.frame())

  min_exons <- combined_min_exons(event_type)

  # Per-transcript stats: exon count and total spliced length
  tx_stats <- dt[, .(
    exon_count = .N,
    spliced    = sum(exon_end - exon_start + 1L)
  ), by = .(transcript_id, gene_id, chr, strand)]

  eligible_tx <- tx_stats[exon_count >= min_exons & spliced >= 500L]

  message(sprintf(
    "  [DIAG] eligible_for_type('%s'): %d / %d transcripts pass filters (>= %d exons, >= 500 bp)",
    event_type, nrow(eligible_tx), nrow(tx_stats), min_exons
  ))

  if (nrow(eligible_tx) == 0L) return(data.frame())

  dt <- dt[transcript_id %in% eligible_tx$transcript_id]

  base_type <- sub("_ES_[12]$", "", event_type)

  if (base_type %in% c("alternative_donor_late",   "alternative_donor_early",
                        "alternative_acceptor_late", "alternative_acceptor_early")) {
    rows <- build_alt_boundary_candidates(dt, base_type, event_type)
  } else if (event_type %in% c("1_ES", "2_ES")) {
    rows <- build_es_candidates(dt, event_type)
  } else if (event_type == "intron_retention") {
    rows <- build_ir_candidates(dt)
  } else {
    message(sprintf("  [DIAG] eligible_for_type: unrecognised event type '%s'", event_type))
    rows <- data.frame()
  }

  message(sprintf(
    "  [DIAG] eligible_for_type('%s'): %d candidate sites found",
    event_type, nrow(rows)
  ))
  rows
}

#' Enforce one event per gene — shuffle first so selection is random.
#' (Caller is responsible for setting the seed before calling this.)
one_event_per_gene <- function(candidates) {
  # FIX Bug 1 (partial): guard against NULL input in addition to zero-row input.
  # The primary fix is upstream (do.call → NULL becomes data.frame()), but this
  # makes the function itself robust to either case.
  if (is.null(candidates) || nrow(candidates) == 0L) return(data.frame())
  candidates <- candidates[sample(nrow(candidates)), ]
  candidates[!duplicated(candidates$gene_id), ]
}

#' Proportional chromosome stratification: sample n rows from candidates,
#' allocating slots to chromosomes in proportion to their share of the pool.
stratify_by_chrom <- function(candidates, n) {
  if (is.null(candidates) || nrow(candidates) == 0L || n == 0L)
    return(if (is.null(candidates)) data.frame() else candidates[integer(0L), ])
  chr_counts <- table(candidates$chr)
  alloc      <- round(n * chr_counts / sum(chr_counts))
  alloc[alloc == 0L] <- 1L
  # Trim back to n if rounding pushed over
  while (sum(alloc) > n) {
    biggest        <- which.max(alloc)
    alloc[biggest] <- alloc[biggest] - 1L
  }
  selected <- lapply(names(alloc), function(ch) {
    pool <- candidates[candidates$chr == ch, ]
    k    <- min(alloc[ch], nrow(pool))
    if (k == 0L) return(pool[integer(0L), ])
    pool[sample(nrow(pool), k), ]
  })
  do.call(rbind, selected)
}

# ── Snakemake glue ────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")

  set.seed(snakemake@params[["seed"]])

  exons <- load_gtf(snakemake@input[["collapsed_gtf"]])

  # Early diagnostic: report what load_gtf returned so problems are immediately visible
  message(sprintf("[INFO] Loaded collapsed GTF: %d features", length(exons)))
  if (!is.null(exons$type)) {
    type_tab <- sort(table(as.character(exons$type)), decreasing = TRUE)
    message("[INFO] Feature type counts: ",
            paste(names(type_tab), type_tab, sep = "=", collapse = ", "))
  } else {
    message("[WARN] No 'type' column found in loaded GTF — assuming all rows are exons")
  }

  # Filter to non-overlapping transcripts to prevent splice-site safety check failures
  message("[INFO] Filtering to non-overlapping transcripts...")
  non_overlapping <- get_non_overlapping_transcripts(exons)
  n_before_filter <- if (!is.null(exons$transcript_id)) {
    length(unique(exons$transcript_id[!is.na(exons$transcript_id)]))
  } else 0L
  exons <- exons[exons$transcript_id %in% non_overlapping]
  n_after_filter <- length(non_overlapping)
  message(sprintf("[INFO] Retained %d / %d transcripts with no exon overlaps",
                  n_after_filter, n_before_filter))

  events    <- snakemake@params[["events"]]
  max_per   <- snakemake@params[["max_per_type"]]
  max_total <- snakemake@params[["max_total"]]

  all_candidates <- list()
  selected_transcripts <- character(0)  # Track transcripts already selected
  for (etype in names(events)) {
    requested <- min(events[[etype]], max_per)
    pool      <- eligible_for_type(exons, etype)
    if (nrow(pool) == 0L) {
      log_warn(sprintf("No eligible candidates for %s — skipping", etype))
      next
    }
    
    # Exclude transcripts already selected for other event types
    pool <- pool[!(pool$transcript_id %in% selected_transcripts), ]
    if (nrow(pool) == 0L) {
      log_warn(sprintf("No eligible candidates for %s after filtering selected transcripts", etype))
      next
    }
    
    selected <- stratify_by_chrom(pool, requested)
    if (nrow(selected) < requested)
      log_warn(sprintf("%s: requested %d, only %d available",
                       etype, requested, nrow(selected)))
    log_info(sprintf("%s: %d / %d candidates selected",
                     etype, nrow(selected), requested))
    
    # Add newly selected transcripts to the tracking set
    selected_transcripts <- c(selected_transcripts, selected$transcript_id)
    
    all_candidates[[etype]] <- selected
  }

  # FIX Bug 1: do.call(rbind, list()) returns NULL, not data.frame().
  # NULL passed to one_event_per_gene causes nrow(NULL) == NULL,
  # making if (NULL == 0L) throw "argument is of length zero".
  candidates <- if (length(all_candidates) == 0L) {
    message("[WARN] No candidates found for any event type — output will be empty")
    data.frame()
  } else {
    do.call(rbind, all_candidates)
  }

  candidates <- one_event_per_gene(candidates)

  if (nrow(candidates) > max_total)
    candidates <- candidates[sample(nrow(candidates), max_total), ]

  log_info(sprintf("Total candidates written: %d", nrow(candidates)))
  write.table(candidates, snakemake@output[[1]],
              sep = "\t", quote = FALSE, row.names = FALSE)

  sink(); sink(type = "message")
  close(log_con)
}