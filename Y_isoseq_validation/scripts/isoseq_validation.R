suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Event type routing constants
# ---------------------------------------------------------------------------

# "up" exon events (dpsi > 0, novel intronic region used as exon in sample):
# validated by checking whether IsoSeq has exonic coverage of the region.
# cAAus: alt acceptor upstream (novel acceptor before canonical)
# cADds: alt donor downstream (novel donor after canonical)
EXON_UP_EVENTS <- c("cAAus",
                     "cADds",
                     "cEI")

# "down" exon events (dpsi < 0, reference exon region absent in sample):
# validated by checking whether the region is intronic in IsoSeq transcripts,
# with an additional splice site coordinate check for the shifted boundary.
# cAAds: alt acceptor downstream (canonical acceptor lost)
# cADus: alt donor upstream (canonical donor lost)
EXON_DOWN_EVENTS <- c("cAAds",
                       "cADus",
                       "cIJ")

# Exon skipping is validated separately: the skipped exon must be present as
# exonic in the control long-read GTF and spanned by an intron in the test GTF.
ES_EVENTS <- c("cES1")

# Intron retention is validated separately regardless of tool.
IR_EVENTS <- c("cIR")

#' Validate novel splicing events against long-read isoform sequencing data
#'
#' @param classified_files Character vector of paths to classified junction files
#' @param isoseq_gtf Path to GTF file from long-read sequencing (e.g., IsoSeq, ONT) - test condition
#' @param output_dir Directory for output files
#' @param tolerance Tolerance in bp for matching splice sites (default: 2)
#' @param isoseq_introns Optional pre-built intronic region index
#' @param control_isoseq_gtf Optional path to control-condition long-read GTF, required for cES1 validation
#' @return Data frame with validation summary statistics
validate_with_longread <- function(classified_files,
                                   isoseq_gtf,
                                   control_isoseq_gtf,
                                   output_dir     = ".",
                                   tolerance      = 2,
                                   isoseq_introns = NULL) {

  cat("=== Long-Read Validation Pipeline ===\n\n")
  cat(sprintf("IsoSeq GTF: %s\n", isoseq_gtf))
  cat(sprintf("Tolerance: %d bp\n", tolerance))
  cat(sprintf("Number of tools to validate: %d\n\n", length(classified_files)))

  # Load IsoSeq annotation
  cat("Loading long-read annotation...\n")
  isoseq_gtf_gr <- rtracklayer::import(isoseq_gtf)
  isoseq_exons  <- isoseq_gtf_gr[isoseq_gtf_gr$type == "exon"]

  cat(sprintf("  Loaded %d exons from %d transcripts\n",
              length(isoseq_exons),
              length(unique(isoseq_exons$transcript_id))))

  # Build IsoSeq junction index
  isoseq_junctions <- build_junction_index_from_gtf(isoseq_exons, tolerance)

  cat(sprintf("  Extracted %d unique junctions\n", nrow(isoseq_junctions)))
  cat(sprintf("  Chromosomes: %s\n\n",
              paste(unique(isoseq_junctions$chr), collapse = ", ")))

  # Build intronic index from test GTF (needed for EXON_DOWN_EVENTS and ES validation)
  if (is.null(isoseq_introns)) {
    isoseq_introns <- build_intronic_index(isoseq_exons)
    cat(sprintf("  Built intronic region index with %d introns from %d transcripts\n",
                nrow(isoseq_introns), length(unique(isoseq_introns$transcript_id))))
  }

  # Load control long-read GTF and build control indices (required)
  cat(sprintf("Loading control long-read annotation: %s\n", control_isoseq_gtf))
  control_gtf_gr         <- rtracklayer::import(control_isoseq_gtf)
  control_isoseq_exons   <- control_gtf_gr[control_gtf_gr$type == "exon"]
  if (length(control_isoseq_exons) == 0) stop("No exons found in control long-read GTF.")
  cat(sprintf("  Loaded %d exons from %d control transcripts\n",
              length(control_isoseq_exons),
              length(unique(control_isoseq_exons$transcript_id))))
  control_junctions      <- build_junction_index_from_gtf(control_isoseq_exons, tolerance)
  control_isoseq_introns <- build_intronic_index(control_isoseq_exons)
  cat(sprintf("  Extracted %d control junctions, %d control introns\n\n",
              nrow(control_junctions), nrow(control_isoseq_introns)))

  # Validate each tool's output
  validation_summary <- list()

  for (i in seq_along(classified_files)) {
    file      <- classified_files[i]
    tool_name <- tools::file_path_sans_ext(basename(file))

    cat(sprintf("=== Validating: %s ===\n", tool_name))

    events    <- fread(file)
    validated <- validate_events(events, isoseq_junctions, isoseq_exons,
                                 tolerance, isoseq_introns,
                                 control_isoseq_exons, control_junctions,
                                 control_isoseq_introns)

    # Add tool name
    validated$tool <- tool_name

    # Calculate summary statistics
    summary_stats <- calculate_validation_stats(validated, tool_name)
    validation_summary[[tool_name]] <- summary_stats

    # Save validated output
    output_file <- file.path(output_dir,
                             paste0(tool_name, "_validated.txt"))
    fwrite(validated, output_file, sep = "\t", quote = FALSE)

    cat(sprintf("  Saved to: %s\n\n", output_file))
  }

  # Combine and save summary
  summary_df  <- bind_rows(validation_summary)
  summary_file <- file.path(output_dir, "validation_summary.txt")
  fwrite(summary_df, summary_file, sep = "\t", quote = FALSE)

  # Print overall summary
  print_validation_summary(summary_df)

  return(summary_df)
}

#' Build intronic region index from exons
#'
#' @param exons GRanges object with exon annotations
#' @return Data table with intronic region coordinates and transcript IDs
build_intronic_index <- function(exons) {

  normalize_chr <- function(x) sub("^chr", "", as.character(x))

  exon_dt <- data.table(
    chr           = normalize_chr(seqnames(exons)),
    start         = start(exons),
    end           = end(exons),
    strand        = as.character(strand(exons)),
    transcript_id = exons$transcript_id
  )

  exon_dt <- exon_dt[order(transcript_id, start)]

  introns_list <- list()

  for (tx_id in unique(exon_dt$transcript_id)) {
    tx_exons <- exon_dt[transcript_id == tx_id]

    if (nrow(tx_exons) < 2) next

    for (j in 1:(nrow(tx_exons) - 1)) {
      exon1        <- tx_exons[j]
      exon2        <- tx_exons[j + 1]
      intron_start <- exon1$end + 1
      intron_end   <- exon2$start - 1

      if (intron_start <= intron_end) {
        introns_list[[length(introns_list) + 1]] <- data.table(
          chr           = exon1$chr,
          intron_start  = intron_start,
          intron_end    = intron_end,
          strand        = exon1$strand,
          transcript_id = tx_id
        )
      }
    }
  }

  if (length(introns_list) == 0) return(data.table())

  introns_dt <- rbindlist(introns_list)
  setkey(introns_dt, chr, strand, intron_start, intron_end)
  return(introns_dt)
}

#' Build junction index from GTF exons
#'
#' @param exons GRanges object with exon annotations
#' @param tolerance Tolerance for matching
#' @return Data table with junction coordinates and transcript IDs
build_junction_index_from_gtf <- function(exons, tolerance) {

  exon_dt <- data.table(
    chr           = sub("^chr", "", as.character(seqnames(exons))),
    start         = start(exons),
    end           = end(exons),
    strand        = as.character(strand(exons)),
    transcript_id = exons$transcript_id
  )

  exon_dt <- exon_dt[order(transcript_id, start)]

  junctions_list <- list()

  for (tx_id in unique(exon_dt$transcript_id)) {
    tx_exons <- exon_dt[transcript_id == tx_id]

    if (nrow(tx_exons) < 2) next

    for (j in 1:(nrow(tx_exons) - 1)) {
      exon1 <- tx_exons[j]
      exon2 <- tx_exons[j + 1]

      junctions_list[[length(junctions_list) + 1]] <- data.table(
        chr           = exon1$chr,
        js            = exon1$end,
        je            = exon2$start,
        strand        = exon1$strand,
        transcript_id = tx_id
      )
    }
  }

  junctions_dt      <- rbindlist(junctions_list)
  junctions_grouped <- junctions_dt[, .(
    transcript_ids = list(unique(transcript_id)),
    n_transcripts  = length(unique(transcript_id))
  ), by = .(chr, js, je, strand)]

  setkey(junctions_grouped, chr, strand, js, je)
  return(junctions_grouped)
}

#' Validate individual events against long-read isoform sequencing data
#'
#' Handles both junction-based tools (columns js/je) and exon-based tools
#' (columns start/end).  Source is detected from column structure:
#'   - Native js/je columns → junction tool
#'   - Native start/end without js/je → exon mode tool
#'
#' Dispatch table (exon mode):
#'   intron_retention               → validate_intron_retention
#'                                      test: exon spans retained region
#'                                      control: canonical junction present
#'   EXON_UP_EVENTS (dpsi > 0)     → validate_exon_coverage
#'                                      test: exon coverage + novel splice site
#'                                      control: canonical splice site or intron
#'   EXON_DOWN_EVENTS (dpsi < 0)   → validate_exon_node_intronic
#'                                      test: intronic containment + novel splice site
#'                                      control: canonical splice site or exon coverage
#'   ES_EVENTS (cES1)              → validate_exon_skipping
#'                                      control: exon spans region + boundaries match junctions; test: intron spans
#'
#' Dispatch table (junction tools):
#'   intron_retention               → validate_intron_retention (as above)
#'   dpsi >= 0 or NA               → validate_junction vs test junctions
#'   dpsi < 0                      → validate_junction vs control junctions
#'
#' @param events Data table of classified events
#' @param isoseq_junctions Long-read junction index (test condition)
#' @param isoseq_exons Long-read exon GRanges (test condition)
#' @param tolerance Tolerance for splice site matching
#' @param isoseq_introns Intronic region index (test condition)
#' @param control_isoseq_exons GRanges of control long-read exons
#' @param control_junctions Long-read junction index (control condition)
#' @param control_isoseq_introns Intronic region index (control condition)
#' @return Data table with isoseq_confirmed and isoseq_transcript_ids columns added
validate_events <- function(events, isoseq_junctions, isoseq_exons,
                            tolerance, isoseq_introns = NULL,
                            control_isoseq_exons   = NULL,
                            control_junctions      = NULL,
                            control_isoseq_introns = NULL) {

  is_exon_mode <- (!"js" %in% colnames(events)) && ("start" %in% colnames(events))

  if (is_exon_mode) {
    cat("  Detected exon-mode input (start/end coordinates).\n")
    events$js <- events$start
    events$je <- events$end
  } else {
    cat("  Detected junction-mode input (js/je coordinates).\n")
  }

  events$isoseq_confirmed      <- FALSE
  events$isoseq_transcript_ids <- NA_character_
  events <- unique(events)

  has_dpsi <- "dpsi" %in% colnames(events)
  has_cluster_id <- "cluster_id" %in% colnames(events)
  has_cEI <- "cEI" %in% colnames(events)

  cat(sprintf("  Total events: %d\n", nrow(events)))
  cat("  Validating events against long-read data...\n")

  for (i in 1:nrow(events)) {
    event      <- events[i, ]
    event_type <- if (!is.na(event$event_type)) event$event_type else "unknown"
    event_desc <- sprintf("%s:%d-%d (%s)", event$chr, event$js, event$je, event$strand)
    dpsi_val   <- if (has_dpsi && !is.na(event$dpsi)) event$dpsi else NA_real_

    if (event_type %in% IR_EVENTS) {
      # Intron retention: same logic for both tool modes.
      # Test: long-read exon must fully span the retained region.
      # Control: canonical splicing junction must be present.
      cat(sprintf("    [%3d/%d] IR     %s\n", i, nrow(events), event_desc))
      result <- validate_intron_retention(
        event$chr, event$js, event$je, event$strand,
        isoseq_exons,
        control_junctions, tolerance
      )

    } else if (is_exon_mode && event_type %in% EXON_UP_EVENTS) {
      cat(sprintf("    [%3d/%d] %-30s %s\n", i, nrow(events), event_type, event_desc))
      result <- validate_exon_coverage(
        event$chr, event$js, event$je, event$strand,
        isoseq_junctions, control_junctions, control_isoseq_introns,
        tolerance, event_type
      )

    } else if (is_exon_mode && event_type %in% EXON_DOWN_EVENTS) {
      cat(sprintf("    [%3d/%d] %-30s %s\n", i, nrow(events), event_type, event_desc))
      if (!is.null(isoseq_introns) && nrow(isoseq_introns) > 0) {
        result <- validate_exon_node_intronic(
          event$chr, event$js, event$je, event$strand,
          isoseq_introns, isoseq_junctions,
          control_junctions, control_isoseq_exons,
          tolerance, event_type
        )
      } else {
        result <- list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_)
      }

    } else if (is_exon_mode && event_type %in% ES_EVENTS) {
      cat(sprintf("    [%3d/%d] %-30s %s\n", i, nrow(events), event_type, event_desc))
      if (!is.null(control_isoseq_exons) && !is.null(isoseq_introns) && nrow(isoseq_introns) > 0) {
        result <- validate_exon_skipping(
          event$chr, event$js, event$je, event$strand,
          control_isoseq_exons, isoseq_introns,
          control_junctions, tolerance
        )
      } else {
        cat(sprintf("      WARNING: cES1 cannot be validated without control GTF; marked NA.\n"))
        result <- list(confirmed = NA, n_transcripts = 0L, transcript_ids = NA_character_)
      }

    } else {
      # Junction-mode tools: dpsi < 0 → losing junction, validate against control.
      #                       dpsi >= 0 or NA → novel junction, validate against test.
      cat(sprintf("    [%3d/%d] SPLICE %s [%s]\n", i, nrow(events), event_desc, event_type))
      result <- validate_junction(
        event$chr, event$js, event$je, event$strand,
        isoseq_junctions, control_junctions,
        tolerance, dpsi_val
      )
    }

    events$isoseq_confirmed[i]      <- result$confirmed
    events$isoseq_transcript_ids[i] <- result$transcript_ids
    if (has_dpsi) events$dpsi[i]    <- event$dpsi
  }

    n_confirmed  <- sum(events$isoseq_confirmed)
  pct_confirmed <- 100 * n_confirmed / nrow(events)

  cat(sprintf("  ============================================\n"))
  cat(sprintf("  Results: %d / %d events confirmed (%.1f%%)\n",
              n_confirmed, nrow(events), pct_confirmed))
  cat(sprintf("  ============================================\n\n"))

  return(events)
}

#' Validate a splice junction against long-read junction data
#'
#' Both boundary coordinates (js and je) must match a long-read junction within
#' the given tolerance.  The junction index queried depends on dpsi direction:
#'   dpsi >= 0 or NA  →  novel/gained junction validated against test junctions
#'   dpsi < 0         →  losing junction validated against control junctions
#'
#' When control_junctions is NULL and dpsi < 0 the function falls back to the
#' test index with a warning rather than returning NA, to preserve backward
#' compatibility when no control GTF has been provided.
#'
#' @param chr Chromosome
#' @param js Junction start coordinate (geometrically lower; donor on +, acceptor on -)
#' @param je Junction end coordinate (geometrically higher; acceptor on +, donor on -)
#' @param strand Strand
#' @param isoseq_junctions Long-read junction index for the test condition
#' @param control_junctions Long-read junction index for the control condition, or NULL
#' @param tolerance Tolerance in bp
#' @param dpsi Delta-PSI value used to select the appropriate index; NA treated as >= 0
#' @return List with confirmed, n_transcripts, transcript_ids
validate_junction <- function(chr, js, je, strand, isoseq_junctions,
                              control_junctions = NULL, tolerance, dpsi = NA_real_) {

  if (length(chr) == 0 || length(js) == 0 || length(je) == 0 || length(strand) == 0 ||
      is.na(chr)  || is.na(js)  || is.na(je)  || is.na(strand)) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  # Select junction index based on dpsi direction.
  # dpsi < 0 → losing junction: must be supported in control.
  # dpsi >= 0 or NA → novel junction: must be supported in test.
  use_control <- !is.na(dpsi) && dpsi < 0
  if (use_control && is.null(control_junctions)) {
    stop("dpsi < 0 event encountered but no control junction index provided. ",
         "Supply control_isoseq_gtf to validate losing junctions.")
  }
  jnc_index <- if (use_control) control_junctions else isoseq_junctions

  chr_clean <- sub("^chr", "", as.character(chr))

  chr_matches <- jnc_index[jnc_index$chr == chr_clean &
                             jnc_index$strand == strand, ]

  if (nrow(chr_matches) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  js_diff          <- abs(chr_matches$js - js)
  je_diff          <- abs(chr_matches$je - je)
  within_tolerance <- (js_diff <= tolerance) & (je_diff <= tolerance)
  matches          <- chr_matches[within_tolerance, ]

  if (nrow(matches) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  all_transcripts <- unique(unlist(matches$transcript_ids))
  all_transcripts <- all_transcripts[!is.na(all_transcripts) & nchar(all_transcripts) > 0]

  if (length(all_transcripts) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  return(list(
    confirmed      = TRUE,
    n_transcripts  = length(all_transcripts),
    transcript_ids = paste(all_transcripts, collapse = ",")
  ))
}

#' Validate a novel exon region (UP events) against long-read data
#'
#' Used for "up" exon events (dpsi > 0): alternative_donor_late,
#' alternative_acceptor_early, cryptic_exon.
#'
#' Three criteria must all be satisfied:
#'
#'   1. Exon coverage in test: at least one IsoSeq exon overlaps the query
#'      region (junction coordinates validated against both boundaries).
#'
#'   2. Novel splice site in test junctions (strand-aware):
#'      Strand logic — on + strand: js (lower coord) = donor, je = acceptor.
#'                      on - strand: js = acceptor, je (higher coord) = donor.
#'      donor_col   = "js" (+) / "je" (-)
#'      acceptor_col = "je" (+) / "js" (-)
#'      five_prime_coord  = js (+) / je (-)   [upstream boundary of contig]
#'      three_prime_coord = je (+) / js (-)   [downstream boundary of contig]
#'
#'      - cADds:    three_prime_coord must match donor_col  in test
#'      - cAAus: five_prime_coord must match acceptor_col in test
#'      - cEI: js → "je" AND je → "js" in test (geometrically invariant)
#'
#'   3. Canonical splice site in control junctions (if control_junctions provided):
#'      - alternative_donor_late:     five_prime_coord must match donor_col  in control
#'      - alternative_acceptor_early: three_prime_coord must match acceptor_col in control
#'      - cryptic_exon:               region must be fully spanned by a control intron
#'   4. Updated terminology (2024-06): "alternative_donor_late" = "cADds", "alternative_acceptor_early" = "cAAus", int
#' @param chr Chromosome
#' @param js Region start (geometrically lower)
#' @param je Region end   (geometrically higher)
#' @param strand Strand
#' @param isoseq_exons IsoSeq exon GRanges (test condition)
#' @param isoseq_junctions Long-read junction index (test condition)
#' @param control_junctions Long-read junction index (control condition), or NULL
#' @param control_isoseq_introns Intronic index from control GTF (for cryptic_exon), or NULL
#' @param tolerance Tolerance in bp for splice site matching
#' @param event_type Event type string
#' @return List with confirmed, n_transcripts, transcript_ids
validate_exon_coverage <- function(chr, js, je, strand,
                                   isoseq_junctions, control_junctions,
                                   control_isoseq_introns,
                                   tolerance, event_type) {

  if (length(chr) == 0 || length(js) == 0 || length(je) == 0 || length(strand) == 0 ||
      is.na(chr)  || is.na(js)  || is.na(je)  || is.na(strand)) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  chr_clean <- sub("^chr", "", as.character(chr))

  # Strand-aware coordinate helpers
  donor_col         <- if (strand == "+") "js" else "je"
  acceptor_col      <- if (strand == "+") "je" else "js"
  five_prime_coord  <- if (strand == "+") js   else je   # upstream contig boundary
  three_prime_coord <- if (strand == "+") je   else js   # downstream contig boundary

  check_site <- function(coord, site_col, jnc) {
    if (nrow(jnc) == 0) return(character(0))
    hits <- jnc[abs(jnc[[site_col]] - coord) <= tolerance, ]
    if (nrow(hits) == 0) return(character(0))
    tx <- unique(unlist(hits$transcript_ids))
    tx[!is.null(tx) & !is.na(tx) & nchar(tx) > 0]
  }

  # ── Step 1: Novel splice site in test junctions ──────────────────────────
  # All UP event types are fully validated by both boundary coordinates:
  # the 80% exon-overlap check is therefore redundant and has been removed.
  jnc_test <- isoseq_junctions[isoseq_junctions$chr == chr_clean &
                                 isoseq_junctions$strand == strand, ]

  if (event_type == "cADds") {
    # Novel donor = three_prime_coord; verify in donor column of test junctions
    site_tx <- check_site(three_prime_coord, donor_col, jnc_test)

  } else if (event_type == "cAAus") {
    # Novel acceptor = five_prime_coord; verify in acceptor column of test junctions
    site_tx <- check_site(five_prime_coord, acceptor_col, jnc_test)

  } else if (event_type == "cEI") {
    # Both splice sites novel; geometrically: js → "je" (upstream junction),
    # je → "js" (downstream junction) — invariant across strands
    tx_5p   <- check_site(js, "je", jnc_test)
    tx_3p   <- check_site(je, "js", jnc_test)
    site_tx <- intersect(tx_5p, tx_3p)

  } else {
    site_tx <- unique(unlist(jnc_test$transcript_ids))
  }

  site_tx <- site_tx[!is.null(site_tx) & !is.na(site_tx) & nchar(site_tx) > 0]

  if (length(site_tx) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  # ── Step 2: Canonical splice site in control junctions ───────────────────
  if (!is.null(control_junctions)) {
    jnc_ctrl <- control_junctions[control_junctions$chr == chr_clean &
                                   control_junctions$strand == strand, ]

    if (event_type == "cAAus") {
      # Canonical acceptor = three_prime_coord; must appear in control acceptor column
      canonical_ok <- length(check_site(three_prime_coord, acceptor_col, jnc_ctrl)) > 0

    } else if (event_type == "cADds") {
      # Canonical donor = five_prime_coord; must appear in control donor column
      canonical_ok <- length(check_site(five_prime_coord, donor_col, jnc_ctrl)) > 0

    } else if (event_type == "cEI") {
      # Cryptic exon should be intronic in control: a control intron must span it
      canonical_ok <- FALSE
      if (!is.null(control_isoseq_introns) && nrow(control_isoseq_introns) > 0) {
        spanning <- control_isoseq_introns[
          control_isoseq_introns$chr    == chr_clean   &
          control_isoseq_introns$strand == strand      &
          control_isoseq_introns$intron_start <= js    &
          control_isoseq_introns$intron_end   >= je
        ]
        canonical_ok <- nrow(spanning) > 0
      }
    } else {
      canonical_ok <- TRUE
    }

    if (!canonical_ok) {
      return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
    }
  }

  return(list(
    confirmed      = TRUE,
    n_transcripts  = length(site_tx),
    transcript_ids = paste(site_tx, collapse = ",")
  ))
}

#' Validate a lost exon region (DOWN events) against long-read data
#'
#' Used for "down" exon events (dpsi < 0): alternative_acceptor_late,
#' alternative_donor_early, exonic_cryptic.
#'
#' Three criteria must all be satisfied (where applicable):
#'
#'   
#'
#'   1. Novel splice site in test junctions (strand-aware):
#'      donor_col   = "js" (+) / "je" (-)
#'      acceptor_col = "je" (+) / "js" (-)
#'      five_prime_coord  = js (+) / je (-)
#'      three_prime_coord = je (+) / js (-)
#'
#'      - alternative_acceptor_late:  three_prime_coord must match acceptor_col in test
#'      - alternative_donor_early:    five_prime_coord  must match donor_col   in test
#'      - exonic_cryptic:             intronic containment only (no unambiguous splice site)
#'
#'   2. Canonical splice site in control junctions (if control_junctions provided):
#'      - alternative_acceptor_late:  five_prime_coord  must match acceptor_col in control
#'      - alternative_donor_early:    three_prime_coord must match donor_col    in control
#'      - exonic_cryptic:             region must have exon coverage in control
#'
#' Transcript IDs reported are from the intersection of containment and site-matched sets.
#'
#' @param chr Chromosome
#' @param js Region start (geometrically lower)
#' @param je Region end   (geometrically higher)
#' @param strand Strand
#' @param isoseq_introns IsoSeq intronic region index (test condition)
#' @param isoseq_junctions Long-read junction index (test condition)
#' @param control_junctions Long-read junction index (control condition), or NULL
#' @param control_isoseq_exons GRanges of control exons (for exonic_cryptic), or NULL
#' @param tolerance Tolerance in bp for splice site matching
#' @param event_type Event type string
#' @return List with confirmed, n_transcripts, transcript_ids
validate_exon_node_intronic <- function(chr, js, je, strand, isoseq_introns,
                                        isoseq_junctions, control_junctions,
                                        control_isoseq_exons,
                                        tolerance, event_type) {

  if (length(chr) == 0 || length(js) == 0 || length(je) == 0 || length(strand) == 0 ||
      is.na(chr)  || is.na(js)  || is.na(je)  || is.na(strand)) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  if (nrow(isoseq_introns) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  normalize_chr <- function(x) sub("^chr", "", as.character(x))
  chr_norm      <- normalize_chr(chr)
  strand_val    <- as.character(strand)

  # Strand-aware coordinate helpers
  donor_col         <- if (strand_val == "+") "js" else "je"
  acceptor_col      <- if (strand_val == "+") "je" else "js"
  five_prime_coord  <- if (strand_val == "+") js   else je
  three_prime_coord <- if (strand_val == "+") je   else js

  check_site <- function(coord, site_col, jnc) {
    if (nrow(jnc) == 0) return(character(0))
    hits <- jnc[abs(jnc[[site_col]] - coord) <= tolerance, ]
    if (nrow(hits) == 0) return(character(0))
    tx <- unique(unlist(hits$transcript_ids))
    tx[!is.null(tx) & !is.na(tx) & nchar(tx) > 0]
  }


  # ── Step 1: Novel splice site in test junctions ───────────────────────────
  jnc_test <- isoseq_junctions[isoseq_junctions$chr    == chr_norm  &
                                 isoseq_junctions$strand == strand_val, ]

  if (event_type == "cAAds") {
    # New late acceptor = three_prime_coord; must appear in acceptor column of test
    site_tx <- check_site(three_prime_coord, acceptor_col, jnc_test)
    confirming_transcripts <- site_tx
    confirming_transcripts <- confirming_transcripts[
    !is.na(confirming_transcripts) &
    nchar(confirming_transcripts) > 0
    ]
  } else if (event_type == "cADus") {
    # New early donor = five_prime_coord; must appear in donor column of test
    site_tx <- check_site(five_prime_coord, donor_col, jnc_test)
    confirming_transcripts <- site_tx
    confirming_transcripts <- confirming_transcripts[
    !is.na(confirming_transcripts) &
    nchar(confirming_transcripts) > 0
    ]
  } else {
    matching_introns <- isoseq_introns[
      isoseq_introns$chr    == chr_norm  &
      isoseq_introns$strand == strand_val &
      intron_start <= js &
      intron_end   >= je
    ]

    if (nrow(matching_introns) == 0) {
      return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
    }

    containment_transcripts <- unique(matching_introns$transcript_id)
    containment_transcripts <- containment_transcripts[!is.na(containment_transcripts)]

    confirming_transcripts <- containment_transcripts

  }

  if (length(confirming_transcripts) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  # ── Step 2: Canonical splice site / exon coverage in control ─────────────
  if (!is.null(control_junctions)) {
    jnc_ctrl <- control_junctions[control_junctions$chr    == chr_norm  &
                                   control_junctions$strand == strand_val, ]

    if (event_type == "cAAds") {
      # Canonical acceptor = five_prime_coord; must appear in control acceptor column
      canonical_ok <- length(check_site(five_prime_coord, acceptor_col, jnc_ctrl)) > 0

    } else if (event_type == "cADus") {
      # Canonical donor = three_prime_coord; must appear in control donor column
      canonical_ok <- length(check_site(three_prime_coord, donor_col, jnc_ctrl)) > 0

    } else if (event_type == "cIJ") {
      # cIJ (intraexonic junction) region should be exonically covered in control
      canonical_ok <- FALSE
      if (!is.null(control_isoseq_exons) && length(control_isoseq_exons) > 0) {
        seqnames_ctrl <- normalize_chr(seqnames(control_isoseq_exons))
        ctrl_ov <- control_isoseq_exons[
          seqnames_ctrl == chr_norm &
          as.character(strand(control_isoseq_exons)) == strand_val &
          start(control_isoseq_exons) <= js &
          end(control_isoseq_exons)   >= je
        ]
        canonical_ok <- length(ctrl_ov) > 0
      }
    } else {
      canonical_ok <- TRUE
    }

    if (!canonical_ok) {
      return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
    }
  }

  return(list(
    confirmed      = TRUE,
    n_transcripts  = length(confirming_transcripts),
    transcript_ids = paste(confirming_transcripts, collapse = ",")
  ))
}

#' Validate exon skipping against control and test long-read data
#'
#' The skipped exon region must satisfy two criteria:
#'   1. Both boundaries match canonical splice sites in the control junctions,
#'      confirming the exon is flanked by annotated junctions (js = exon end, je = exon start).
#'   2. It is fully spanned by a single intron in the test long-read data,
#'      confirming the exon is skipped in the test condition.
#'
#' @param chr Chromosome
#' @param js Skipped exon start
#' @param je Skipped exon end
#' @param strand Strand
#' @param control_isoseq_exons GRanges of exons from the control long-read GTF (not used directly; only for function signature compatibility)
#' @param isoseq_introns Intronic region index built from the test long-read GTF
#' @param control_junctions Long-read junction index (control condition), or NULL
#' @param tolerance Tolerance in bp for splice site matching
#' @return List with confirmed, n_transcripts, transcript_ids
validate_exon_skipping <- function(chr, js, je, strand,
                                   control_isoseq_exons, isoseq_introns,
                                   control_junctions = NULL, tolerance = 2) {

  if (length(chr) == 0 || length(js) == 0 || length(je) == 0 || length(strand) == 0 ||
      is.na(chr)  || is.na(js)  || is.na(je)  || is.na(strand)) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  chr_clean  <- sub("^chr", "", as.character(chr))
  strand_val <- as.character(strand)

  # ── Step 1: Boundary matching in control junctions ────────────────────────
  # Both boundaries must match canonical splice sites in control:
  # js = exon end (donor on +, acceptor on -), je = exon start (acceptor on +, donor on -)
  boundary_ok <- TRUE
  if (!is.null(control_junctions)) {
    jnc_ctrl <- control_junctions[control_junctions$chr    == chr_clean  &
                                   control_junctions$strand == strand_val, ]

    # Check if js matches an exon end (donor position) and je matches an exon start (acceptor position)
    js_matches <- jnc_ctrl[abs(jnc_ctrl$js - js) <= tolerance, ]
    je_matches <- jnc_ctrl[abs(jnc_ctrl$je - je) <= tolerance, ]

    boundary_ok <- (nrow(js_matches) > 0) && (nrow(je_matches) > 0)

    if (!boundary_ok) {
      return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
    }
  }

  # ── Step 2: Exon fully spanned by an intron in test ───────────────────────
  if (nrow(isoseq_introns) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  spanning_introns <- isoseq_introns[
    isoseq_introns$chr == chr_clean &
    isoseq_introns$strand == strand_val &
    isoseq_introns$intron_start <= js &
    isoseq_introns$intron_end   >= je
  ]

  if (nrow(spanning_introns) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  confirming_transcripts <- unique(spanning_introns$transcript_id)
  confirming_transcripts <- confirming_transcripts[!is.na(confirming_transcripts)]

  if (length(confirming_transcripts) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  return(list(
    confirmed      = TRUE,
    n_transcripts  = length(confirming_transcripts),
    transcript_ids = paste(confirming_transcripts, collapse = ",")
  ))
}

#' Validate intron retention against long-read data
#'
#' Two criteria must both be satisfied:
#'
#'   1. Test condition: a long-read exon must completely span [js, je],
#'      confirming that the intron is treated as exonic in the test condition.
#'
#'   2. Control condition (if control_junctions provided): a long-read junction
#'      with js ≈ js_event and je ≈ je_event must exist in the control data,
#'      confirming that the intron is normally spliced in the control.
#'      This check is strand-agnostic because js < je holds for all junctions
#'      regardless of strand, and a canonical junction always has its lower
#'      exon end as js and upper exon start as je.
#'
#' @param chr Chromosome
#' @param js Retained intron start (geometrically lower)
#' @param je Retained intron end (geometrically higher)
#' @param strand Strand
#' @param isoseq_exons IsoSeq exon GRanges (test condition)
#' @param control_junctions Long-read junction index (control condition), or NULL
#' @param tolerance Tolerance in bp for control junction matching
#' @return List with confirmed, n_transcripts, transcript_ids
validate_intron_retention <- function(chr, js, je, strand,
                                      isoseq_exons,
                                      control_junctions, tolerance) {

  if (length(chr) == 0 || length(js) == 0 || length(je) == 0 || length(strand) == 0 ||
      is.na(chr)  || is.na(js)  || is.na(je)  || is.na(strand)) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  normalize_chr  <- function(x) sub("^chr", "", as.character(x))
  chr_clean      <- normalize_chr(chr)
  seqnames_clean <- normalize_chr(seqnames(isoseq_exons))

  # Test condition IsoSeq exon must fully span the retained intron
  spanning <- isoseq_exons[
    seqnames_clean == chr_clean &
    as.character(strand(isoseq_exons)) == strand &
    start(isoseq_exons) <= js &
    end(isoseq_exons)   >= je
  ]

  if (length(spanning) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  confirming_transcripts <- unique(as.character(spanning$transcript_id))
  confirming_transcripts <- confirming_transcripts[!is.na(confirming_transcripts)]

  if (length(confirming_transcripts) == 0) {
    return(list(confirmed = FALSE, n_transcripts = 0L, transcript_ids = NA_character_))
  }

  return(list(
    confirmed      = TRUE,
    n_transcripts  = length(confirming_transcripts),
    transcript_ids = paste(confirming_transcripts, collapse = ",")
  ))
}

#' Calculate validation statistics
#'
#' @param validated Data table with validation results
#' @param tool_name Name of the tool
#' @return Data frame with summary statistics
calculate_validation_stats <- function(validated, tool_name) {

  overall <- data.frame(
    tool          = tool_name,
    event_type    = "ALL",
    total         = nrow(validated),
    confirmed     = sum(validated$isoseq_confirmed),
    pct_confirmed = 100 * sum(validated$isoseq_confirmed) / nrow(validated),
    stringsAsFactors = FALSE
  )

  by_event <- validated %>%
    group_by(event_type) %>%
    summarize(
      total         = n(),
      confirmed     = sum(isoseq_confirmed),
      pct_confirmed = 100 * sum(isoseq_confirmed) / n(),
      .groups       = "drop"
    ) %>%
    mutate(tool = tool_name) %>%
    select(tool, event_type, total, confirmed, pct_confirmed)

  return(bind_rows(overall, by_event))
}

#' Print validation summary
#'
#' @param summary_df Summary data frame
print_validation_summary <- function(summary_df) {

  cat("\n=== Validation Summary ===\n\n")

  overall <- summary_df %>%
    filter(event_type == "ALL") %>%
    arrange(desc(pct_confirmed))

  cat("Overall confirmation rates by tool:\n")
  cat(sprintf("%-20s %10s %10s %10s\n", "Tool", "Total", "Confirmed", "% Confirmed"))
  cat(strrep("-", 55), "\n")

  for (i in 1:nrow(overall)) {
    cat(sprintf("%-20s %10d %10d %9.1f%%\n",
                overall$tool[i],
                overall$total[i],
                overall$confirmed[i],
                overall$pct_confirmed[i]))
  }

  cat("\n")

  by_event <- summary_df %>%
    filter(event_type != "ALL") %>%
    group_by(event_type) %>%
    summarize(
      total         = sum(total),
      confirmed     = sum(confirmed),
      pct_confirmed = 100 * sum(confirmed) / sum(total)
    ) %>%
    arrange(desc(pct_confirmed))

  cat("Confirmation rates by event type (all tools):\n")
  cat(sprintf("%-35s %10s %10s %10s\n", "Event Type", "Total", "Confirmed", "% Confirmed"))
  cat(strrep("-", 70), "\n")

  for (i in 1:nrow(by_event)) {
    cat(sprintf("%-35s %10d %10d %9.1f%%\n",
                by_event$event_type[i],
                by_event$total[i],
                by_event$confirmed[i],
                by_event$pct_confirmed[i]))
  }

  cat("\n")
}

#' Generate validation report plots
#'
#' @param summary_df Summary data frame
#' @param output_dir Output directory
generate_validation_plots <- function(summary_df, output_dir = ".") {

  library(ggplot2)

  # Plot 1: Overall confirmation by tool
  overall <- summary_df %>% filter(event_type == "ALL")

  p1 <- ggplot(overall, aes(x = reorder(tool, pct_confirmed),
                             y = pct_confirmed, fill = tool)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", pct_confirmed)),
              hjust = -0.2, size = 3) +
    coord_flip() +
    ylim(0, 100) +
    theme_minimal() +
    labs(title = "IsoSeq Confirmation Rate by Tool",
         x = "Tool", y = "% Confirmed") +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "validation_by_tool.pdf"),
         p1, width = 8, height = 6)

  # Plot 2: Confirmation by event type
  by_event <- summary_df %>%
    filter(event_type != "ALL") %>%
    group_by(event_type) %>%
    summarize(
      total         = sum(total),
      confirmed     = sum(confirmed),
      pct_confirmed = 100 * sum(confirmed) / sum(total)
    )

  p2 <- ggplot(by_event, aes(x = reorder(event_type, pct_confirmed),
                              y = pct_confirmed)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", pct_confirmed)),
              hjust = -0.2, size = 3) +
    coord_flip() +
    ylim(0, 100) +
    theme_minimal() +
    labs(title = "IsoSeq Confirmation Rate by Event Type",
         x = "Event Type", y = "% Confirmed")

  ggsave(file.path(output_dir, "validation_by_event.pdf"),
         p2, width = 10, height = 6)

  # Plot 3: Heatmap of tool vs event type
  heatmap_data <- summary_df %>%
    filter(event_type != "ALL") %>%
    select(tool, event_type, pct_confirmed)

  p3 <- ggplot(heatmap_data, aes(x = event_type, y = tool,
                                  fill = pct_confirmed)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.0f%%", pct_confirmed)), size = 3) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                         midpoint = 50, limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "IsoSeq Confirmation: Tool vs Event Type",
         x = "Event Type", y = "Tool", fill = "% Confirmed")

  ggsave(file.path(output_dir, "validation_heatmap.pdf"),
         p3, width = 12, height = 6)

  cat(sprintf("Plots saved to %s\n", output_dir))
}

# ============================================================================
# Main Snakemake Execution Block
# ============================================================================

tryCatch({

  # Get input/output from Snakemake
  input_file  <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]

  # Get parameters
  isoseq_gtf <- snakemake@params[["isoseq_gtf"]]
  tolerance  <- snakemake@params[["tolerance"]]
  if (is.null(tolerance) || length(tolerance) == 0) {
    tolerance <- 2
  } else {
    tolerance <- as.integer(tolerance)
  }
  if (length(tolerance) == 0 || is.na(tolerance)) {
    stop("tolerance parameter is missing or invalid")
  }

  output_dir <- snakemake@params[["output_dir"]]
  if (is.null(output_dir) || length(output_dir) == 0) {
    output_dir <- dirname(snakemake@output[[1]])
  }

  # Optional log file
  log_file <- NULL
  if (!is.null(snakemake@log) && length(snakemake@log) > 0) {
    log_file <- snakemake@log[[1]]
    log_con  <- file(log_file, open = "wt")
    sink(log_con, split = FALSE)
    sink(log_con, type = "message")
  }

  cat("=== IsoSeq Validation ===\n\n")
  cat(sprintf("Input file: %s\n", input_file))
  cat(sprintf("IsoSeq GTF: %s\n", isoseq_gtf))
  cat(sprintf("Output file: %s\n", output_file))
  cat(sprintf("Tolerance: %d bp\n\n", tolerance))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load IsoSeq annotation
  cat("Loading long-read annotation...\n")
  isoseq_gtf_gr <- rtracklayer::import(isoseq_gtf)
  isoseq_exons  <- isoseq_gtf_gr[isoseq_gtf_gr$type == "exon"]
  if (length(isoseq_exons) == 0) stop("No exons found in IsoSeq GTF")

  cat(sprintf("  Loaded %d exons from %d transcripts\n",
              length(isoseq_exons),
              length(unique(isoseq_exons$transcript_id))))

  # Build IsoSeq junction index
  isoseq_junctions <- build_junction_index_from_gtf(isoseq_exons, tolerance)

  cat(sprintf("  Extracted %d unique junctions\n", nrow(isoseq_junctions)))
  cat(sprintf("  Chromosomes: %s\n\n",
              paste(unique(isoseq_junctions$chr), collapse = ", ")))

  # Build intronic index from test GTF (needed for EXON_DOWN_EVENTS and ES validation)
  isoseq_introns <- build_intronic_index(isoseq_exons)
  cat(sprintf("  Built intronic region index with %d introns from %d transcripts\n",
              nrow(isoseq_introns), length(unique(isoseq_introns$transcript_id))))

  # Control long-read GTF — required for two-sided validation
  control_isoseq_gtf <- snakemake@params[["control_isoseq_gtf"]]
  if (is.null(control_isoseq_gtf) || nchar(control_isoseq_gtf) == 0)
    stop("control_isoseq_gtf must be supplied in snakemake@params for two-sided validation.")

  cat(sprintf("\nLoading control long-read annotation: %s\n", control_isoseq_gtf))
  control_gtf_gr         <- rtracklayer::import(control_isoseq_gtf)
  control_isoseq_exons   <- control_gtf_gr[control_gtf_gr$type == "exon"]
  if (length(control_isoseq_exons) == 0) stop("No exons found in control long-read GTF.")
  cat(sprintf("  Loaded %d exons from %d control transcripts\n",
              length(control_isoseq_exons),
              length(unique(control_isoseq_exons$transcript_id))))
  control_junctions      <- build_junction_index_from_gtf(control_isoseq_exons, tolerance)
  control_isoseq_introns <- build_intronic_index(control_isoseq_exons)
  cat(sprintf("  Extracted %d control junctions, %d control introns\n",
              nrow(control_junctions), nrow(control_isoseq_introns)))

  # Load and validate classified events
  cat(sprintf("\nLoading classified events: %s\n", basename(input_file)))
  events <- fread(input_file)
  cat(sprintf("  Loaded %d events\n", nrow(events)))

  validated <- validate_events(events, isoseq_junctions, isoseq_exons,
                               tolerance, isoseq_introns,
                               control_isoseq_exons, control_junctions,
                               control_isoseq_introns)

  # Confirmation summary
  cat("\n=== Confirmation Summary ===\n")
  cat(sprintf("Data type of isoseq_confirmed: %s\n", class(validated$isoseq_confirmed)))
  n_true  <- sum(validated$isoseq_confirmed == TRUE,  na.rm = TRUE)
  n_false <- sum(validated$isoseq_confirmed == FALSE, na.rm = TRUE)
  n_na    <- sum(is.na(validated$isoseq_confirmed))
  cat(sprintf("TRUE  : %d\n", n_true))
  cat(sprintf("FALSE : %d\n", n_false))
  cat(sprintf("NA    : %d\n", n_na))
  cat(sprintf("Total : %d\n", nrow(validated)))

  # Save validated output
  cat(sprintf("\nWriting validated output to: %s\n", output_file))
  fwrite(validated, output_file, sep = "\t", quote = FALSE)

  # Verify write
  cat("\n=== Verification of Written File ===\n")
  written_data <- fread(output_file)
  cat(sprintf("Rows written: %d\n", nrow(written_data)))
  cat(sprintf("Columns: %s\n", paste(colnames(written_data), collapse = ", ")))

  if ("isoseq_confirmed" %in% colnames(written_data)) {
    written_true <- sum(written_data$isoseq_confirmed == TRUE, na.rm = TRUE)
    cat(sprintf("Confirmed (written): %d / %d (%.1f%%)\n",
                written_true, nrow(written_data),
                100 * written_true / nrow(written_data)))

    summary_file <- file.path(output_dir, "validation_counts.txt")
    cat(sprintf("Tool: %s\n",              basename(input_file)),   file = summary_file, append = FALSE)
    cat(sprintf("Total events: %d\n",      nrow(written_data)),      file = summary_file, append = TRUE)
    cat(sprintf("Confirmed events: %d\n",  written_true),            file = summary_file, append = TRUE)
    cat(sprintf("Confirmation rate: %.1f%%\n",
                100 * written_true / nrow(written_data)),            file = summary_file, append = TRUE)
    cat(sprintf("Output file MD5: %s\n",
                system(sprintf("md5sum %s | awk '{print $1}'", output_file), intern = TRUE)),
        file = summary_file, append = TRUE)
  }
  cat("================================\n\n")

  cat("✓ Validation complete!\n")

  if (!is.null(log_file)) {
    sink()
    sink(type = "message")
    close(log_con)
  }

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  if (exists("log_file") && !is.null(log_file)) {
    sink()
    sink(type = "message")
    close(log_con)
  }
  quit(status = 1)
})