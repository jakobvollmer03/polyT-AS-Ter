#!/usr/bin/env Rscript
## This script classifies a set of novel junctions (either previously identified by get_junctions or 
## provided by a splice detection tool). The classifier assumes that everything that it gets is confirmed to be 
## novel. Therefore it classifies:
## - All junctions whose input event_type column equals "intron_retention" are classified as INTRON RETENTION.
##   They pass through the boundary-matching logic so that gene_id, transcript_id and reference coordinates
##   are populated, but their event_type is fixed to "intron_retention" regardless of what the classifier
##   would otherwise assign.
## - All junctions for which BOTH boundaries match an annotated exon end and exon start of a CANONICAL
##   junction are DISCARDED (not emitted in the output).
## - All junctions for which the genomically FIRST boundary MATCHES, but not the second one, as ALTERNATIVE ACCEPTOR
## - All junctions for which the SECOND boundary MATCHES, but not the first one, as ALTERNATIVE DONORS
## - All junctions for which BOTH boundaries MATCH exon boundaries, but NOT of ADJACENT exons, as EXON SKIPPING (ES)
## - All junctions for which NEITHER boundaries MATCH, as COMPLEX
## - Additionally, the classifier will add the amount of exons skipped, if any, to every event type. It will also 
## add information on whether non-matching boundaries are genomically upstream (early) or downstream (late) of 
## the respective canonical boundary
##
## Input table format: chr, intron_start, intron_end, strand, event_type  [+ optional extra columns]
## event_type column: junctions with event_type == "intron_retention" run through the classifier normally
##                    to obtain coordinates/gene info, but always exit with event_type = "intron_retention";
##                    all other values (including NA) are classified without any override.

## Example: alternative_acceptor_late_ES_2 (first boundary matches a canonical exon end, 2 exons are skipped and the
## second boundary is within the range of the third canonical exon)


# Splice Junction Event Type Classifier - REFACTORED VERSION
# Pending refinements
# - add option for dpsi to seperate IR and ES events
# - 

# Required libraries quietly
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(data.table)
  library(stringr)
}) 

# Get input/output from Snakemake
# Support both named and positional parameter access for backward compatibility
if (!is.null(snakemake@input[["junction_file"]])) {
  junction_file <- snakemake@input[["junction_file"]]
} else if (!is.null(snakemake@input[["novel_features"]])) {
  junction_file <- snakemake@input[["novel_features"]]
} else {
  junction_file <- snakemake@params[[1]]  # Fallback to positional
}

gtf_file <- if (!is.null(snakemake@params[["gtf_file"]])) {
  snakemake@params[["gtf_file"]]
} else if (!is.null(snakemake@params[["reference_gtf"]])) {
  snakemake@params[["reference_gtf"]]
} else {
  snakemake@params[[2]]  # Fallback to positional
}

output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]
tolerance <- if (!is.null(snakemake@params[["tolerance"]])) {
  as.numeric(snakemake@params[["tolerance"]])
} else {
  2  
}
isoseq <- isTRUE(snakemake@params[["isoseq"]])

# Redirect output to log file
if (!is.null(log_file)) {
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = FALSE)
  sink(log_con, type = "message")
}

load_reference <- function(gtf_file) {
  cat("Loading reference annotation...\n")
  
  # Import GTF
  gtf <- rtracklayer::import(gtf_file)
  
  # Filter for exons only
  exons <- gtf[gtf$type == "exon"]
  
  cat(sprintf("Loaded %d exons from %d transcripts\n", 
              length(exons), 
              length(unique(exons$transcript_id))))
  
  return(exons)
}

#' Load junction table
load_junctions <- function(junction_file, isoseq = FALSE) {
  cat("Loading junction table...\n")
  
  # Read file, skipping comment lines
  junctions <- as.data.table(read.table(junction_file, 
                                  header = TRUE, 
                                  comment.char = "#",
                                  fill = TRUE, 
                                  stringsAsFactors = FALSE,
                                  row.names = NULL))
  
  if (!isoseq) {

  # Format: chr, intron_start, intron_end, strand, [read_count]
  # Convert to exon boundaries: js = intron_start - 1, je = intron_end + 1

  # Handle case where column names might already be set
    if (ncol(junctions) >= 4 && 
        !all(colnames(junctions)[1:4] == c("chr", "intron_start", "intron_end", "strand"))) {
      
      colnames(junctions)[1:4] <- c("chr", "intron_start", "intron_end", "strand")
      
      junctions$intron_start <- as.integer(junctions$intron_start)
      junctions$intron_end <- as.integer(junctions$intron_end)
      
      junctions <- junctions %>%
        mutate(
          js = as.integer(intron_start - 1),  # last base of upstream exon
          je = as.integer(intron_end + 1)     # first base of downstream exon
        )
      
    } else if (ncol(junctions) >= 4) {
      
      junctions <- junctions %>%
        mutate(
          chr = V1,
          intron_start = as.integer(V2),
          intron_end = as.integer(V3),
          strand = V4,
          js = as.integer(V2 - 1),
          je = as.integer(V3 + 1)
        )
      
    } else {
      stop("Junction file must have at least 4 columns: chr, intron_start, intron_end, strand")
    }

  } else if (isoseq) {

    # In isoseq mode input tables have varying formats, but always contain columns: chr, start, end, strand.
    # start/end are in intron coordinates, convert to exon boundaries
    if (ncol(junctions) >= 4) {
      
      junctions <- junctions %>%
        mutate(
          chr = as.character(chr),
          strand = as.character(strand),
          js = as.integer(start),  
          je = as.integer(end)
        )
      
      # Handle dpsi column if it exists
      if ("dpsi" %in% colnames(junctions)) {
        junctions <- junctions %>%
          mutate(dpsi = as.numeric(dpsi))
      }
      
    } else {
      stop("IsoSeq format requires at least 4 columns: chr, strand, start, end")
    }
  }

  # Ensure required columns exist
  if (!all(c("chr", "js", "je", "strand") %in% colnames(junctions))) {
    stop("Missing required columns after format conversion")
  }

  # Handle event_type column (5th input column).
  # When the input file carries a proper header, read.table() will have already
  # named the column "event_type".  In the header-less branch the column will be
  # present as "V5".  In either case we normalise to "event_type" here.
  if ("feature_type" %in% colnames(junctions)) {
    setnames(junctions, "feature_type", "event_type")
  } 

  if (!"event_type" %in% colnames(junctions) && !"feature_type" %in% colnames(junctions)) {
    if ("V5" %in% colnames(junctions)) {
      setnames(junctions, "V5", "event_type")
    } else if (ncol(junctions) >= 5) {
      # Rename whatever is in position 5 as a last resort
      setnames(junctions, colnames(junctions)[5], "event_type")
    } else {
      junctions$event_type <- NA_character_
    }
  }
  junctions$event_type <- as.character(junctions$event_type)
  if ("event_type" %in% colnames(junctions)) {
    cat(sprintf("Event types found in input:\n"))
    print(table(junctions$event_type, useNA = "ifany"))
  } else {
    cat("  No event types found in input. Input columns:\n")
    print(colnames(junctions))
  }
  
  # Validate coordinates
  if (any(junctions$js >= junctions$je, na.rm = TRUE)) {
    warning("Some junctions have js >= je. This indicates a data problem.")
    cat("Examples of problematic junctions:\n")
    print(head(junctions[junctions$js >= junctions$je, ]))
  }
  
  # Add unique junction ID
  junctions$junction_id <- paste(junctions$chr, 
                                 junctions$js, 
                                 junctions$je, 
                                 junctions$strand, 
                                 sep = "_")
  
  cat(sprintf("Loaded %d junctions\n", nrow(junctions)))
  cat(sprintf("  Chromosomes: %s\n", paste(unique(junctions$chr), collapse = ", ")))
  cat(sprintf("  Strands: %s\n", paste(unique(junctions$strand), collapse = ", ")))
  
  return(junctions)
}

#' Build exon boundary index for fast lookup
#'
#' @param exons GRanges object with exon annotations
#' @return List with exon start and end indices
build_exon_index <- function(exons) {
  cat("Building exon boundary index...\n")
  
  # Create data table for fast lookups
  dt <- data.table(
    chr = as.character(seqnames(exons)),
    strand = as.character(strand(exons)),
    exon_start = start(exons),
    exon_end = end(exons),
    gene_id = exons$gene_id,
    transcript_id = exons$transcript_id,
    exon_number = exons$exon_number
  )
  
  # Index by exon ends (for js matching)
  setkey(dt, chr, strand, exon_end)
  end_index <- copy(dt)
  
  # Index by exon starts (for je matching)
  setkey(dt, chr, strand, exon_start)
  start_index <- copy(dt)
  
  return(list(exon_end = end_index, exon_start = start_index))
}

#' Find exons matching a specific boundary position
#'
#' @param chr Chromosome
#' @param pos Position
#' @param strand Strand
#' @param index Exon boundary index
#' @param boundary_type "exon_end" or "exon_start"
#' @param tolerance Maximum distance to consider as a match (default: 0, exact match only)
#' @return Data table of matching exons
find_matching_exons <- function(chr, pos, strand, index, boundary_type = "exon_end", tolerance = 0) {
  
  idx <- if (boundary_type == "exon_end") index$exon_end else index$exon_start
  
  tryCatch({
    # Get the boundary column name (exon_start or exon_end)
    boundary_col <- if (boundary_type == "exon_end") "exon_end" else "exon_start"
    
    # Get all records for this chr and strand using regular filtering
    # (Can't use binary search with partial keys)
    all_records <- idx[chr == chr & strand == strand]
    
    if (is.null(all_records) || nrow(all_records) == 0) {
      return(data.table())
    }
    
    # Calculate distances to find matches within tolerance
    all_records[, distance := abs(get(boundary_col) - pos)]
    matches <- all_records[distance <= tolerance]
    
    if (nrow(matches) == 0) {
      return(data.table())
    }
    
    # Filter out NA results and remove temporary distance column
    matches <- matches[!is.na(gene_id)]
    matches[, distance := NULL]
    
    if (nrow(matches) == 0) {
      return(data.table())
    }
    
    return(matches)
  }, error = function(e) {
    # Return empty data table on error
    return(data.table())
  })
}

#' Classify a single junction
#'
#' @param junction Single row from junction table
#' @param index Exon boundary index
#' @param exons_gr GRanges object with all exons
#' @param tolerance Maximum distance for boundary matching
#' @return Data frame with classification results
classify_junction <- function(junction, index, exons_gr, tolerance = 2) {
  
  chr <- junction$chr
  js <- junction$js      # junction start (geometrically lower coordinate)
  je <- junction$je      # junction end (geometrically higher coordinate)
  strand <- junction$strand
  
  # Check if this is a pre-labeled intron retention event
  is_prelabeled_ir <- !is.na(junction$event_type) && junction$event_type == "intron_retention"
  is_prelabeled_en <- !is.na(junction$event_type) && junction$event_type == "exon_node"


  # Initialize result (needed by every return path)
  result <- data.frame(
    junction_id = junction$junction_id,
    chr = chr,
    js = js,
    je = je,
    strand = strand,
    event_type = "novel",
    gene_id = NA_character_,
    transcript_id = NA_character_,
    ref_js = NA_integer_,
    ref_je = NA_integer_,
    confidence = "0",
    dpsi = if ("dpsi" %in% colnames(junction)) junction$dpsi else NA_real_,
    stringsAsFactors = FALSE
  )

  # Find exact matches for js and je (with tolerance)
  js_matches <- find_matching_exons(chr, js, strand, index, "exon_end", tolerance)
  je_matches <- find_matching_exons(chr, je, strand, index, "exon_start", tolerance)
  
  js_matches$exon_number <- as.integer(js_matches$exon_number)
  je_matches$exon_number <- as.integer(je_matches$exon_number)
  
  # ------------------------------------------
  # FORCE handling of pre-labeled IR junctions
  # ------------------------------------------
  if (is_prelabeled_ir) {

    # Try boundary matches first
    if (nrow(js_matches) > 0) {
      result$gene_id <- js_matches$gene_id[1]
      result$transcript_id <- js_matches$transcript_id[1]
    }

    if (nrow(je_matches) > 0 && is.na(result$gene_id)) {
      result$gene_id <- je_matches$gene_id[1]
      result$transcript_id <- je_matches$transcript_id[1]
    }

    # Try overlapping exon (covers partial IR cases)
    if (is.na(result$gene_id)) {
      overlapping_exons <- exons_gr[
        seqnames(exons_gr) == chr &
        strand(exons_gr) == strand &
        start(exons_gr) <= je &
        end(exons_gr) >= js
      ]

      if (length(overlapping_exons) > 0) {
        result$gene_id <- overlapping_exons$gene_id[1]
        result$transcript_id <- overlapping_exons$transcript_id[1]
        result$ref_js <- start(overlapping_exons[1])
        result$ref_je <- end(overlapping_exons[1])
      }
    }

    # Fallback: nearest gene
    if (is.na(result$gene_id)) {
      all_exons_chr <- exons_gr[
        seqnames(exons_gr) == chr &
        strand(exons_gr) == strand
      ]

      if (length(all_exons_chr) > 0) {
        nearest_idx <- which.min(abs(start(all_exons_chr) - js))
        result$gene_id <- all_exons_chr$gene_id[nearest_idx]
        result$transcript_id <- all_exons_chr$transcript_id[nearest_idx]
      }
    }

    # Always set coordinates (even if imperfect)
    result$ref_js <- js
    result$ref_je <- je

    # Check dpsi: if negative, the sequence is spliced out in test group → exonic_cryptic
    if (!is.na(result$dpsi) && result$dpsi < 0) {
      result$event_type <- "exonic_cryptic"
      result$confidence <- "pre-labeled_IR_dpsi<0"
    } else {
      result$event_type <- "intron_retention"
      result$confidence <- "pre-labeled"
    }

    return(result)
  }





  # Special handling for pre-labeled intron retention events:
  # If this is marked as IR in the input, it will always be output as IR
  # after boundary matching to populate gene/transcript info
  if (is_prelabeled_ir) {
    # Attempt to find matching genes/transcripts
    if (nrow(js_matches) > 0 && nrow(je_matches) > 0) {
      # Both boundaries match - use first matching gene
      common_genes <- intersect(js_matches$gene_id, je_matches$gene_id)
      if (length(common_genes) > 0) {
        result$gene_id <- common_genes[1]
        result$transcript_id <- js_matches$transcript_id[1]
        result$ref_js <- js
        result$ref_je <- je
        if (!is.na(result$dpsi) && result$dpsi < 0) {
          result$event_type <- "exonic_cryptic"
          result$confidence <- "pre-labeled_IR_dpsi<0"
        } else {
          result$event_type <- "intron_retention"
          result$confidence <- "pre-labeled"
        }
        return(result)
      }
    }
    
    if (nrow(js_matches) > 0) {
      # js matches - use this gene
      result$gene_id <- js_matches$gene_id[1]
      result$transcript_id <- js_matches$transcript_id[1]
      result$ref_js <- js
      result$ref_je <- je
      if (!is.na(result$dpsi) && result$dpsi < 0) {
        result$event_type <- "exonic_cryptic"
        result$confidence <- "pre-labeled_IR_dpsi<0"
      } else {
        result$event_type <- "intron_retention"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    if (nrow(je_matches) > 0) {
      # je matches - use this gene
      result$gene_id <- je_matches$gene_id[1]
      result$transcript_id <- je_matches$transcript_id[1]
      result$ref_js <- js
      result$ref_je <- je
      if (!is.na(result$dpsi) && result$dpsi < 0) {
        result$event_type <- "exonic_cryptic"
        result$confidence <- "pre-labeled_IR_dpsi<0"
      } else {
        result$event_type <- "intron_retention"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    # No boundary matches - try to find overlapping exons
    overlapping_exons <- exons_gr[seqnames(exons_gr) == chr &
                                   strand(exons_gr) == strand &
                                   start(exons_gr) <= js &
                                   end(exons_gr) >= je]
    
    if (length(overlapping_exons) > 0) {
      # Junction is within an exon
      result$gene_id <- overlapping_exons$gene_id[1]
      result$transcript_id <- overlapping_exons$transcript_id[1]
      result$ref_js <- start(overlapping_exons[1])
      result$ref_je <- end(overlapping_exons[1])
      if (!is.na(result$dpsi) && result$dpsi < 0) {
        result$event_type <- "exonic_cryptic"
        result$confidence <- "pre-labeled_IR_dpsi<0"
      } else {
        result$event_type <- "intron_retention"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    # Last resort: try to find nearest gene
    all_exons_chr <- exons_gr[seqnames(exons_gr) == chr & strand(exons_gr) == strand]
    if (length(all_exons_chr) > 0) {
      nearest_idx <- which.min(abs(start(all_exons_chr) - js))
      result$gene_id <- all_exons_chr$gene_id[nearest_idx]
      result$transcript_id <- all_exons_chr$transcript_id[nearest_idx]
      result$ref_js <- js
      result$ref_je <- je
      if (!is.na(result$dpsi) && result$dpsi < 0) {
        result$event_type <- "exonic_cryptic"
        result$confidence <- "pre-labeled_IR_dpsi<0"
      } else {
        result$event_type <- "intron_retention"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    # Even if no exons found, still mark as IR or exonic_cryptic based on dpsi
    if (!is.na(result$dpsi) && result$dpsi < 0) {
      result$event_type <- "exonic_cryptic"
      result$confidence <- "pre-labeled_IR_dpsi<0"
    } else {
      result$event_type <- "intron_retention"
      result$confidence <- "pre-labeled"
    }
    return(result)
  }
  
  # Case 1: Both js and je match exon boundaries
  if (nrow(js_matches) > 0 && nrow(je_matches) > 0) {
    
    # Check if they're from the same transcript
    common_transcripts <- intersect(js_matches$transcript_id, 
                                   je_matches$transcript_id)
    
    if (length(common_transcripts) > 0) {
      # Check each common transcript to see if exons are adjacent
      is_annotated <- FALSE
      for (tx_id in common_transcripts) {
        # Get exons for this transcript
        tx_js_exons <- js_matches[js_matches$transcript_id == tx_id, ]
        tx_je_exons <- je_matches[je_matches$transcript_id == tx_id, ]
        
        # Check if exon numbers are present and valid
        if (all(!is.na(tx_js_exons$exon_number)) && 
            all(!is.na(tx_je_exons$exon_number))) {
          
          # Check if any combination has adjacent exons
          for (js_enum in tx_js_exons$exon_number) {
            for (je_enum in tx_je_exons$exon_number) {
              if (je_enum - js_enum == 1) {
                is_annotated <- TRUE
                # Both boundaries match adjacent exons - this is a canonical junction
                # (it would have been caught in pre-labeled IR handling above if marked as IR)
                result$event_type    <- "canonical_discarded"
                result$confidence    <- "1"
                result$transcript_id <- tx_id
                result$gene_id       <- tx_js_exons$gene_id[1]
                result$ref_js        <- js
                result$ref_je        <- je
                # overwrite event type and confidence if this was pre-labeled as exon node
                if (is_prelabeled_en) {
                  result$event_type <- "exon_node"
                  result$confidence <- "pre-labeled"
                }
                return(result)
              }
            }
          }
        }
      }
    }
    
    # Check for exon skipping: js and je from same gene but non-adjacent exons
    common_genes <- intersect(js_matches$gene_id, je_matches$gene_id)
    
    if (length(common_genes) > 0) {
      gene_id <- common_genes[1]
      
      # Get all exons for this gene on this chromosome and strand
      gene_exons <- exons_gr[exons_gr$gene_id == gene_id & 
                             seqnames(exons_gr) == chr &
                             strand(exons_gr) == strand]
      
      if (length(gene_exons) == 0) {
        result$event_type <- "novel_junction"
        result$confidence <- "1"
        return(result)
      }
      
      # Count skipped exons (those completely within js-je interval)
      # An exon is skipped if: exon_start > js AND exon_end < je
      skipped <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
      # get ref_je by finding start of first downstream exon of js
      downstream_exons <- gene_exons[start(gene_exons) > js]
      if (length(downstream_exons) > 0) {
        ref_je <- min(start(downstream_exons))
      } else {
        ref_je <- NA_integer_
      }
      # since js might not match an annotated exon end, get the closest exon end to js for ref_js
      nearby_exon <- gene_exons[which.min(abs(end(gene_exons) - js))]
      ref_js <- end(nearby_exon)

      if (length(skipped) > 0) {
        n_skipped <- length(unique(skipped))  # Count unique exons
        
        result$event_type <- ifelse(n_skipped == 1, 
                                   "1_ES", 
                                   paste0(n_skipped, "_ES")) # MES and IR of certain isoforms cannot be distinguished purely from junction data
        result$gene_id <- gene_id
        result$transcript_id <- js_matches$transcript_id[1]
        result$confidence <- "1"
        result$ref_js <- ref_js
        result$ref_je <- ref_je
        # overwrite event type and confidence if this was pre-labeled as exon node
        if (is_prelabeled_en) {
          result$event_type <- "exon_node"
          result$confidence <- "pre-labeled"
        }
        return(result)
      } else {
        # Both boundaries match annotated exon boundaries of the same gene with no
        # exons skipped. This is a canonical junction and should be discarded.
        result$event_type    <- "canonical_discarded"
        result$confidence    <- "1"
        result$gene_id       <- gene_id
        result$transcript_id <- js_matches$transcript_id[1]
        result$ref_js        <- ref_js
        result$ref_je        <- ref_je
        # overwrite event type and confidence if this was pre-labeled as exon node
        if (is_prelabeled_en) {
          result$event_type <- "exon_node"
          result$confidence <- "pre-labeled"
        }
        return(result)
      }
    }
  }
  
  # Case 2: js matches, je doesn't (alternative je)
  if (nrow(js_matches) > 0 && nrow(je_matches) == 0) {
    
    gene_id <- js_matches$gene_id[1]
    transcript_id <- js_matches$transcript_id[1]
    
    # Get all exons for this gene on this chromosome and strand
    gene_exons <- exons_gr[exons_gr$gene_id == gene_id & 
                           seqnames(exons_gr) == chr &
                           strand(exons_gr) == strand]
    
    # Find the canonical je for this js
    canon_je <- find_canonical_partner(chr, js, strand, exons_gr, gene_id, 
                                      partner_type = "start")
    # Check if an exon has both boundaries between js and je - if so, it is both an alternative acceptor/donor and exon skipping, which should be added to the event type
    overlapping_exons <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
    if (length(overlapping_exons) > 0) {
      ES <- paste0("_ES_", length(unique(overlapping_exons)))
    } else {
      ES <- ""
    }


    if (!is.na(canon_je)) {
      if (length(overlapping_exons) == 0) {
        # Determine position relative to canonical
        position_plus <- ifelse(je < canon_je, "early", "late")
        position_minus <- ifelse(je > canon_je, "early", "late")
      }
      # if exons are skipped, position_je must be determined relative to the exon downstream of the last on that was skipped.
      if (length(overlapping_exons) > 0) {
        # Find the downstream exon of the last skipped exon
        last_skipped <- max(end(overlapping_exons))
        downstream_exon <- gene_exons[start(gene_exons) > last_skipped]
        if (length(downstream_exon) > 0) {
          ref_je <- start(downstream_exon[1])
          position_plus <- ifelse(je < ref_je, "early", "late")
          position_minus <- ifelse(je > ref_je, "early", "late")
        } else {
          ref_je <- NA_integer_
        }
      } else {
        ref_je <- canon_je
      }
      # Biological interpretation based on strand
      if (strand == "+") {
        # Plus strand: je is acceptor
        result$event_type <- paste0("alternative_acceptor_", position_plus, ES)
      } else {
        # Minus strand: je is donor
        result$event_type <- paste0("alternative_donor_", position_minus, ES)
      }
      # Since js might not exactly match an annotated exon end, get the closest exon end to js for ref_js
      nearby_exon <- gene_exons[which.min(abs(end(gene_exons) - js))]
      ref_js <- end(nearby_exon)

      result$gene_id <- gene_id
      result$transcript_id <- transcript_id
      result$ref_js <- ref_js
      if (length(overlapping_exons) > 0) {
        result$ref_je <- ref_je
      } else {
        result$ref_je <- canon_je
      }
      result$confidence <- "1"
      # overwrite event type and confidence if this was pre-labeled as exon node
      if (is_prelabeled_en) {
        result$event_type <- "exon_node"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    # Fallback: check if je is near a gene exon
    result <- classify_alternative_boundary(chr, js, je, strand, gene_id, 
                                           exons_gr, result, "je")
    if (result$event_type != "novel") {
      return(result)
    }
  }
  
  # Case 3: je matches, js doesn't (alternative js)
  if (nrow(js_matches) == 0 && nrow(je_matches) > 0) {
    
    gene_id <- je_matches$gene_id[1]
    transcript_id <- je_matches$transcript_id[1]
    
    # Get all exons for this gene on this chromosome and strand
    gene_exons <- exons_gr[exons_gr$gene_id == gene_id & 
                           seqnames(exons_gr) == chr &
                           strand(exons_gr) == strand]
    
    # Find the canonical js for this je
    canon_js <- find_canonical_partner(chr, je, strand, exons_gr, gene_id,
                                      partner_type = "end")
    # Check if an exon has both boundaries between js and je - if so, it is both an alternative acceptor/donor and exon skipping, which should be added to the event type
    overlapping_exons <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
    if (length(overlapping_exons) > 0) {
      ES <- paste0("_ES_", length(unique(overlapping_exons)))
      # If exons are skipped, position_js must be determined relative to the exon upstream of the first on that was skipped.
      first_skipped <- min(start(overlapping_exons))
      upstream_exon <- gene_exons[end(gene_exons) < first_skipped]
      if (length(upstream_exon) > 0) {
        ref_js <- end(upstream_exon[length(upstream_exon)])
      } else {
        ref_js <- NA_integer_
      }
    } else {
      ES <- ""
    }

    if (!is.na(canon_js)) {
      if (length(overlapping_exons) == 0) {
      # Determine position relative to canonical
      position_plus <- ifelse(js < canon_js, "early", "late")
      position_minus <- ifelse(js > canon_js, "early", "late")
      }
      # Determine position relative to canonical based on skipped exons if present
    if (length(overlapping_exons) > 0) {
      position_plus <- ifelse(js < ref_js, "early", "late")
      position_minus <- ifelse(js > ref_js, "early", "late")
    }
      # Biological interpretation based on strand
      if (strand == "+") {
        # Plus strand: js is donor
        result$event_type <- paste0("alternative_donor_", position_plus, ES)
      } else {
        # Minus strand: js is acceptor
        result$event_type <- paste0("alternative_acceptor_", position_minus, ES)
      }
      
      # Since je might not exactly match an annotated exon start, get the closest exon start to je for ref_je
      nearby_exon <- gene_exons[which.min(abs(start(gene_exons) - je))]
      ref_je <- start(nearby_exon)
      
      result$gene_id <- gene_id
      result$transcript_id <- transcript_id
      if (length(overlapping_exons) > 0) {
        result$ref_js <- ref_js
      } else {
         result$ref_js <- canon_js
      }
      result$ref_je <- ref_je
      result$confidence <- "1"
      # overwrite event type and confidence if this was pre-labeled as exon node
      if (is_prelabeled_en) {
        result$event_type <- "exon_node"
        result$confidence <- "pre-labeled"
      }
      return(result)
    }
    
    # Fallback: check if js is near a gene exon
    result <- classify_alternative_boundary(chr, js, je, strand, gene_id,
                                           exons_gr, result, "js")
    if (result$event_type != "novel") {
      return(result)
    }
  }
  
  # Case 4: Neither matches - check for intron retention or novel junction
  result <- classify_novel_junction(junction, exons_gr, result, tolerance)
  
  # FINAL CHECK: If this was pre-labeled as exon_node, force it back to exon_node
  # This ensures that all pre-labeled exon_node events are output as exon_node
  # regardless of the classification logic above
  if (is_prelabeled_en) {
    result$event_type <- "exon_node"
    result$confidence <- "pre-labeled"
  }
  
  return(result)
}

#' Find canonical partner boundary for a given exon boundary
#'
#' @param chr Chromosome
#' @param pos Position of known boundary
#' @param strand Strand
#' @param exons_gr GRanges with exons
#' @param gene_id Gene ID to search within
#' @param partner_type "start" (looking for exon start) or "end" (looking for exon end)
#' @return Position of canonical partner or NA
find_canonical_partner <- function(chr, pos, strand, exons_gr, gene_id,
                                  partner_type = "start") {
  
  gene_exons <- exons_gr[exons_gr$gene_id == gene_id & 
                         seqnames(exons_gr) == chr &
                         strand(exons_gr) == strand]
  
  if (length(gene_exons) == 0) return(NA)
  
  # Sort by position
  gene_exons <- gene_exons[order(start(gene_exons))]
  
  if (partner_type == "start") {
    # We have an exon end (js), need to find next exon's start (je)
    next_exons <- gene_exons[start(gene_exons) > pos]
    if (length(next_exons) > 0) {
      return(start(next_exons[1]))
    }
  } else {
    # We have an exon start (je), need to find previous exon's end (js)
    prev_exons <- gene_exons[end(gene_exons) < pos]
    if (length(prev_exons) > 0) {
      return(end(prev_exons[length(prev_exons)]))
    }
  }
  
  return(NA)
}

#' Classify alternative boundary when near but not exact match
#'
#' @param chr Chromosome
#' @param js Junction start
#' @param je Junction end
#' @param strand Strand
#' @param gene_id Gene ID
#' @param exons_gr All exons
#' @param result Current result dataframe
#' @param alt_boundary Which boundary is alternative ("js" or "je")
#' @return Updated result
classify_alternative_boundary <- function(chr, js, je, strand, gene_id,
                                         exons_gr, result, alt_boundary) {
  
  is_prelabeled_en <- !is.na(result$event_type) && result$event_type == "exon_node"
  gene_exons <- exons_gr[exons_gr$gene_id == gene_id &
                         seqnames(exons_gr) == chr &
                         strand(exons_gr) == strand]
  
  if (length(gene_exons) == 0) return(result)
  
  if (alt_boundary == "je") {
    # js matches, je is alternative - find nearest exon start
    nearby_exon <- gene_exons[which.min(abs(start(gene_exons) - je))]
    ref_je <- start(nearby_exon)
    
    # Also need to find the canonical js (closest exon end to js)
    nearby_js_exon <- gene_exons[which.min(abs(end(gene_exons) - js))]
    ref_js <- end(nearby_js_exon)

    # Check if there's an exon completely between js and je - if so, it is both an alternative acceptor/donor and exon skipping, which should be added to the event type
    overlapping_exons <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
    if (length(overlapping_exons) > 0) {
      ES <- paste0("_ES_", length(unique(overlapping_exons)))
    } else {      
      ES <- ""
    }

    
    if (abs(je - ref_je) > tolerance) {
      position_plus <- ifelse(je < ref_je, "early", "late")
      position_minus <- ifelse(je > ref_je, "early", "late")
      
      # Biological interpretation
      if (strand == "+") {
        result$event_type <- paste0("alternative_acceptor_", position_plus, ES)
      } else {
        result$event_type <- paste0("alternative_donor_", position_minus, ES)
      }
      
      result$gene_id <- gene_id
      result$ref_js <- ref_js
      result$ref_je <- ref_je
      result$confidence <- "2"
      # overwrite event type and confidence if this was pre-labeled as exon node
      if (is_prelabeled_en) {
        result$event_type <- "exon_node"
        result$confidence <- "pre-labeled"
      }
    }
    
  } else {
    # je matches, js is alternative - find nearest exon end
    nearby_exon <- gene_exons[which.min(abs(end(gene_exons) - js))]
    ref_js <- end(nearby_exon)
    
    # Check if there's an exon completely between js and je - if so, it is both an alternative acceptor/donor and exon skipping, which should be added to the event type
    overlapping_exons <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
    if (length(overlapping_exons) > 0) {
      ES <- "_ES"
    } else {      
      ES <- ""
    }

    # Only classify if reasonably close (within tolerance)
    if (abs(js - ref_js) < tolerance) {
      position_plus <- ifelse(js < ref_js, "early", "late")
      position_minus <- ifelse(js > ref_js, "early", "late")
      
      # Biological interpretation
      if (strand == "+") {
        result$event_type <- paste0("alternative_donor_", position_plus, ES)
      } else {
        result$event_type <- paste0("alternative_acceptor_", position_minus, ES)
      }
      
      result$gene_id <- gene_id
      result$transcript_id <- nearby_exon$transcript_id[1]
      result$ref_js <- ref_js
      result$ref_je <- ref_je
      result$confidence <- "2"
      # overwrite event type and confidence if this was pre-labeled as exon node
      if (is_prelabeled_en) {
        result$event_type <- "exon_node"
        result$confidence <- "pre-labeled"
      }
    }
  }
  
  return(result)
}

#' Classify novel junctions (neither boundary matches annotation)
#'
#' @param junction Junction data
#' @param exons_gr GRanges with exons
#' @param result Initial result data frame
#' @param tolerance Maximum distance for boundary matching
#' @return Updated result data frame
classify_novel_junction <- function(junction, exons_gr, result, tolerance = 2) {
  
  chr <- junction$chr
  js <- junction$js
  je <- junction$je
  strand <- junction$strand
  
  # Check if this was pre-labeled as an exon node
  is_prelabeled_en <- !is.na(result$event_type) && result$event_type == "exon_node"
  
  # Initialize canon_js and canon_je to NA (they may be set later in conditional blocks)
  canon_js <- NA_integer_
  canon_je <- NA_integer_
  
  # Check if junction falls within an annotated exon (potential cryptic intron)
  overlapping_exons <- exons_gr[seqnames(exons_gr) == chr &
                                 strand(exons_gr) == strand &
                                 start(exons_gr) <= js &
                                 end(exons_gr) >= je]

  # get corresponding transcript ID and ref_js/ref_je for the overlapping exon if it exists
  
  if (length(overlapping_exons) > 0) {
    result$event_type <- "exonic_cryptic"
    result$gene_id <- overlapping_exons$gene_id[1]
    result$transcript_id <- overlapping_exons$transcript_id[1]
    result$ref_js <- start(overlapping_exons[1])
    result$ref_je <- end(overlapping_exons[1])
    result$confidence <- "2"
    # overwrite event type and confidence if this was pre-labeled as exon node
    if (is_prelabeled_en) {
      result$event_type <- "exon_node"
      result$confidence <- "pre-labeled"
    }
    return(result)
  }
  
  # Check if junction connects two exons from the same gene (complex includes any junctions
  # that share neither donor nor acceptor with the reference but still connect exons of the same gene)
  gene_exons <- exons_gr[seqnames(exons_gr) == chr &
                         strand(exons_gr) == strand]
  
  downstream_exons <- gene_exons[start(gene_exons) >= js]
  upstream_exons <- gene_exons[end(gene_exons) <= je]
  
  if (length(upstream_exons) > 0 && length(downstream_exons) > 0) {
    common_genes <- intersect(upstream_exons$gene_id, downstream_exons$gene_id)
    
    if (length(common_genes) > 0) {
      canon_js <- find_canonical_partner(chr, js, strand, exons_gr, common_genes[1], "end")
      canon_je <- find_canonical_partner(chr, je, strand, exons_gr, common_genes[1], "start")
      
      # Check if coordinates are within tolerance of canonical
      js_within_tolerance <- !is.na(canon_js) && abs(js - canon_js) <= tolerance
      je_within_tolerance <- !is.na(canon_je) && abs(je - canon_je) <= tolerance
      # If either or both coordinates are within tolerance of canonical, don't classify as complex
      if (js_within_tolerance || je_within_tolerance) {
        if (js_within_tolerance && je_within_tolerance) {
          # Both boundaries are within tolerance of a canonical junction.
          # Pre-labeled IR events are confirmed here; all others are discarded.
          result$event_type <- "canonical_discarded"
          result$confidence <- "2"
          result$gene_id <- common_genes[1]
          result$ref_js <- canon_js
          result$ref_je <- canon_je
        } else if (js_within_tolerance && !je_within_tolerance) {
          # js within tolerance, je is alternative
          result$event_type <- ifelse(je < canon_je, "alternative_acceptor_early", "alternative_acceptor_late")
          result$gene_id <- common_genes[1]
          result$ref_js <- canon_js
          result$ref_je <- canon_je
          result$confidence <- "2"
        } else if (!js_within_tolerance && je_within_tolerance) {
          # je within tolerance, js is alternative
          result$event_type <- ifelse(js < canon_js, "alternative_donor_early", "alternative_donor_late")
          result$gene_id <- common_genes[1]
          result$ref_js <- canon_js
          result$ref_je <- canon_je
          result$confidence <- "2"
        }
        return(result)
      }
      
      # Find the actual closest canonical boundaries for position classification
      nearby_exon <- gene_exons[which.min(abs(start(gene_exons) - je))]
      ref_je <- start(nearby_exon)
    
      # Also need to find the canonical js (closest exon end to js)
      nearby_js_exon <- gene_exons[which.min(abs(end(gene_exons) - js))]
      ref_js <- end(nearby_js_exon)
      
      if (!is.na(ref_js)) {
        # Determine position relative to canonical
        position_plus <- ifelse(js < ref_js, "early", "late")
        position_minus <- ifelse(js > ref_js, "early", "late")
        position_js <- ifelse(strand == "+", position_plus, position_minus)
      } else {
        position_js <- NA_character_
      }
      
      if (!is.na(ref_je)) {
        # Determine position relative to canonical
        position_plus <- ifelse(je < ref_je, "early", "late")
        position_minus <- ifelse(je > ref_je, "early", "late")
        position_je <- ifelse(strand == "+", position_plus, position_minus)
      } else {
        position_je <- NA_character_
      }

      # Check if any exons are completely between js and je - if so, it is both a complex event and exon skipping, which should be added to the event type
      overlapping_exons <- gene_exons[start(gene_exons) > js & end(gene_exons) < je]
      if (length(overlapping_exons) > 0) {
        ES <- paste0("_ES_", length(unique(overlapping_exons)))
      } else {
        ES <- ""
      }

      if (length(common_genes) > 0) {
        gene_id <- common_genes[1]
        ifelse(strand == "+",
               result$event_type <- paste0("complex_", position_js, "_", position_je, ES),
               result$event_type <- paste0("complex_", position_je, "_", position_js, ES))
        result$gene_id <- gene_id
        result$transcript_id <- gene_exons$transcript_id[1]
        result$confidence <- "2"
        result$ref_js <- ref_js
        result$ref_je <- ref_je
        if (is_prelabeled_en) {
          result$event_type <- "exon_node"
          result$confidence <- "pre-labeled"
        }

        # Safely extract ref_js and ref_je
        #upstream_gene_exons <- upstream_exons[upstream_exons$gene_id == gene_id]
        #downstream_gene_exons <- downstream_exons[downstream_exons$gene_id == gene_id]
        
        #if (length(upstream_gene_exons) > 0) {
        #  result$ref_js <- end(upstream_gene_exons[1])
        #} else {
        #  result$ref_js <- end(downstream_gene_exons[1])
        #}
        #if (length(downstream_gene_exons) > 0) {
        #  result$ref_je <- start(downstream_gene_exons[1])
        #} else {
        #  result$ref_je <- start(upstream_gene_exons[1])
        #}
        
        return(result)
      }
    }
  }
  
  # Otherwise, completely novel
  result$event_type <- "novel_junction"
  result$confidence <- "3"
  result$ref_js <- canon_js
  result$ref_je <- canon_je
  result$gene_id <- NA_character_
  result$transcript_id <- NA_character_
  
  return(result)
}

#' Main classification function
#'
#' @param junction_file Path to junction file
#' @param gtf_file Path to GTF annotation
#' @param output_file Output file path
#' @param isoseq Logical indicating if the input is in IsoSeq format
#' @param tolerance Maximum distance for boundary matching
#' @param cores Number of cores for parallel processing
classify_all_junctions <- function(junction_file, gtf_file, output_file,
                                  isoseq = isoseq, tolerance = 2, cores = 1) {
  
  cat("=== Splice Junction Event Type Classifier ===\n\n")
  cat(sprintf("Tolerance for coordinate matching: %d bp\n\n", tolerance))
  
  # Load data
  exons_gr <- load_reference(gtf_file)
  junctions <- load_junctions(junction_file, isoseq = isoseq)
  
  # Build index
  index <- build_exon_index(exons_gr)
  
  # Classify junctions
  cat("\nClassifying junctions...\n")
  total_prelabeled_ir <- sum(junctions$event_type == "intron_retention", na.rm = TRUE)
  cat(sprintf("Pre-labeled intron retention events: %d\n", total_prelabeled_ir))
  total_prelabeled_es <- sum(junctions$event_type == "exon_node", na.rm = TRUE)
  cat(sprintf("Pre-labeled exon skipping events: %d\n", total_prelabeled_es))
  
  if (cores > 1) {
    library(parallel)
    cl <- makeCluster(cores)
    
    # Export functions AND data
    clusterExport(cl, c("classify_junction", "find_matching_exons", 
                      "find_canonical_partner", "classify_alternative_boundary",
                      "classify_novel_junction", "index", "exons_gr", "tolerance"),
                  envir = environment())
    
    # Load required libraries on each worker
    clusterEvalQ(cl, {
      library(GenomicRanges)
      library(dplyr)
      library(data.table)
    })
    
    results <- pbapply::pblapply(1:nrow(junctions), function(i) {
      classify_junction(junctions[i, ], index, exons_gr, tolerance)
    }, cl = cl)
    
    stopCluster(cl)
  } else {
    results <- list()
    pb <- txtProgressBar(min = 0, max = nrow(junctions), style = 3)
    for (i in 1:nrow(junctions)) {
      results[[i]] <- classify_junction(junctions[i, ], index, exons_gr, tolerance)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  # Combine results
  results_df <- bind_rows(results) #%>% 
    #select(junction_id, chr, js, je, strand, event_type, gene_id, transcript_id,
     #      ref_js, ref_je)
  
  # Remove junctions that matched canonical splice sites (both boundaries annotated).
  # These were previously mis-classified as intron_retention; they are now discarded.
  n_canonical <- sum(results_df$event_type == "canonical_discarded", na.rm = TRUE)
  if (n_canonical > 0) {
    cat(sprintf("\nDiscarding %d canonical junction(s) (both boundaries match annotated exon boundaries).\n",
                n_canonical))
    results_df <- results_df[!(results_df$event_type == "canonical_discarded" &results_df$confidence != "pre-labeled"), ]
  }
  n_novel <- sum(results_df$event_type == "novel_junction", na.rm = TRUE)
  if (n_novel > 0) {
    warning(sprintf("\n%d novel junction(s) detected.\n", n_novel))
  }

  # Summary statistics
  cat("\n\n=== Classification Summary ===\n")
  print(table(results_df$event_type))
  
  # Write output
  cat(sprintf("\nWriting results to %s\n", output_file))
  fwrite(results_df, output_file, sep = "\t", quote = FALSE)
  
  cat("\nDone!\n")
  
  return(results_df)
}

# Main execution with Snakemake
tryCatch({
  cat("Starting junction classification...\n")
  cat(sprintf("Input junction file: %s\n", junction_file))
  cat(sprintf("GTF file: %s\n", gtf_file))
  cat(sprintf("Output file: %s\n", output_file))
  cat("\n")
  
  # Run classification (adjust cores based on threads if available)
  cores <- if (exists("snakemake") && !is.null(snakemake@threads)) {
    as.numeric(snakemake@threads)
  } else {
    1
  }
  
  results <- classify_all_junctions(
    junction_file = junction_file,
    gtf_file = gtf_file,
    output_file = output_file,
    isoseq = isoseq,
    tolerance = tolerance,
    cores = cores
  )
  
  cat("\n✓ Classification complete!\n")
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Close log file
if (!is.null(log_file)) {
  sink()
  sink(type = "message")
  close(log_con)
}