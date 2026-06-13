# function for bidirectional filtering between a strintie gtf and a set of differentially spliced exon contigs 
# Used to compare accuracy current SpliCeAT pipeline against single DS detection tools on the transcript level 


whippet_in     <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r1_ds_detection/whippet/whippet_delta_psi.diff"
stringtie_gtf  <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/stringtie_assemblies/stringtie_assembly.gtf"
reference_gtf  <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r0_get_ref/Mus_musculus/Mus_musculus_GRCm39_115_chr.gtf"
output_dir     <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/masterlists"
novel_gtf      <- "/mnt/gtklab01/dbsv771/SpliCeAT_res/CTX_mm_1205/Meg_version/results/r2_augment_transcriptome__event_mode/whippet_stringtie_novel_transcripts.gtf"
min_aad_length <- 5  # minimum AA/AD node width (bp) — below this, boundary matching is unreliable
boundary_tol   <- 2   # splice site coordinate tolerance (bp)

suppressMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(data.table)
    library(dplyr)
    library(stringr)
    library(tidyr)
})

# Only CE, AA, AD, RI are retained — these are the event types shared across
# junction-based tools (MAJIQ, LeafCutter) and Whippet. AF, AL, TS, TE are
# excluded: AF/AL require transcript-boundary evidence; TS/TE are
# transcription-initiation events invisible to junction-counting tools.
cat("WHIPPET: Reading and filtering differential PSI file...\n")

whippet_raw <- readLines(whippet_in)
whippet_raw <- whippet_raw[whippet_raw != "" & !grepl("^#", whippet_raw)]
whippet_df  <- as.data.frame(do.call(rbind, strsplit(whippet_raw, "\t", fixed = TRUE)))
colnames(whippet_df) <- c("Gene", "Node", "Coord", "Strand", "Type",
                           "Psi_A", "Psi_B", "DeltaPsi", "Probability",
                           "Complexity", "Entropy")

whippet_events <- whippet_df %>%
    filter(grepl("^-?[0-9.]+", DeltaPsi) & grepl("^[0-9.]+", Probability)) %>%
    mutate(across(c(DeltaPsi, Probability), as.numeric)) %>%
    filter(abs(DeltaPsi) >= 0.2 & Probability >= 0.95 &
               Type %in% c("CE","AA","AD","RI")) %>%
    separate(Coord, into = c("chr", "coords"), sep = ":") %>%
    separate(coords, into = c("start", "end"), sep = "-") %>%
    mutate(
        start        = as.numeric(start),
        end          = as.numeric(end),
        genes        = str_remove(str_remove(Gene, "^gene:"), "\\.[0-9]+$"),
        genes        = ifelse(genes == "" | genes == ".", NA_character_, genes),
        tool         = "Whippet",
        feature_type = ifelse(Type == "RI", "intron_retention", Type),
        lsv_id       = paste0(chr, ":", start, "-", end, ":", Strand)
    ) %>%
    dplyr::select(gene_id = genes, chr, strand = Strand, start, end,
                  tool, feature_type, lsv_id, dpsi = DeltaPsi, prob = Probability)

cat(sprintf("WHIPPET: %d events pass filters (CE/AA/AD/RI, |dPSI|>=0.2, prob>=0.95).\n",
            nrow(whippet_events)))

outfile1 <- file.path(output_dir, "whippet_lsvs_w_event_types.tsv")
fwrite(whippet_events, outfile1, sep = "\t")

cat("WHIPPET: Building GRanges and splitting by event type and dPSI direction...\n")

whippet_gr <- GRanges(
    seqnames     = whippet_events$chr,
    ranges       = IRanges(start = whippet_events$start, end = whippet_events$end),
    strand       = whippet_events$strand,
    feature_type = whippet_events$feature_type,
    dpsi         = whippet_events$dpsi,
    lsv_id       = whippet_events$lsv_id
)

ce_pos <- whippet_gr[whippet_gr$feature_type == "CE"               & whippet_gr$dpsi > 0]
ce_neg <- whippet_gr[whippet_gr$feature_type == "CE"               & whippet_gr$dpsi < 0]
ri_pos <- whippet_gr[whippet_gr$feature_type == "intron_retention"  & whippet_gr$dpsi > 0]
ri_neg <- whippet_gr[whippet_gr$feature_type == "intron_retention"  & whippet_gr$dpsi < 0]
aa_pos <- whippet_gr[whippet_gr$feature_type == "AA"               & whippet_gr$dpsi > 0]
aa_neg <- whippet_gr[whippet_gr$feature_type == "AA"               & whippet_gr$dpsi < 0]
ad_pos <- whippet_gr[whippet_gr$feature_type == "AD"               & whippet_gr$dpsi > 0]
ad_neg <- whippet_gr[whippet_gr$feature_type == "AD"               & whippet_gr$dpsi < 0]

# Extract exons directly from the imported GRanges — no TxDb needed.
# TxDb construction (makeTxDbFromGRanges) is the dominant runtime bottleneck
# and provides no benefit here: we only need the flat exon set for coverage,
# and reference introns derived from consecutive exon pairs for section 5b.
cat("TRUTH SET: Processing reference GTF...\n")

official_ref      <- rtracklayer::import(reference_gtf)
ref_exons_flat <- official_ref[
    official_ref$type == "exon" &
    as.character(strand(official_ref)) %in% c("+", "-") &
    #if biotypes are available, restrict to protein_coding — otherwise keep all (e.g. for mouse)
    if (!is.null(official_ref$transcript_biotype)) {
        official_ref$transcript_biotype == "protein_coding"
    } else {
        TRUE
    }
]
ref_exons_reduced <- reduce(ref_exons_flat, ignore.strand = FALSE)
cat(sprintf("TRUTH SET: %d reference exon ranges reduced to %d non-overlapping blocks.\n",
            length(ref_exons_flat), length(ref_exons_reduced)))

# Reference introns derived from consecutive exon pairs per transcript.
# These are needed in section 5b to pre-clear MSTRG introns that exactly match
# a reference intron — such introns represent canonical splicing and must not
# be subjected to the novel-intronic-contig filter. Without this pre-clearance,
# reduce() merging of alternative reference isoforms causes canonical MSTRG
# introns to appear to span alternative exons from other isoforms, falsely
# flagging virtually all transcripts as unsupported in multi-isoform genes.
ref_exon_dt <- data.table(
    chr           = as.character(seqnames(ref_exons_flat)),
    start         = start(ref_exons_flat),
    end           = end(ref_exons_flat),
    strand        = as.character(strand(ref_exons_flat)),
    transcript_id = ref_exons_flat$transcript_id
)[order(transcript_id, start)]

ref_intron_dt <- ref_exon_dt[, {
    n <- .N
    if (n < 2L) data.table(chr=character(), intron_start=integer(),
                            intron_end=integer(), strand=character())
    else data.table(
        chr          = chr[seq_len(n - 1L)],
        intron_start = end[seq_len(n - 1L)] + 1L,
        intron_end   = start[2L:n] - 1L,
        strand       = strand[seq_len(n - 1L)]
    )
}, by = transcript_id][intron_start <= intron_end]

ref_introns_gr <- GRanges(
    seqnames = ref_intron_dt$chr,
    ranges   = IRanges(ref_intron_dt$intron_start, ref_intron_dt$intron_end),
    strand   = ref_intron_dt$strand
)
cat(sprintf("TRUTH SET: %d reference introns derived from exon structure.\n",
            length(ref_introns_gr)))

# Extract exons directly; derive introns from consecutive exon pairs per
# transcript (replicating GenomicFeatures::intronsByTranscript without TxDb).
cat("STRINGTIE: Parsing GTF and deriving introns from exon structure...\n")

gtf_raw   <- rtracklayer::import(stringtie_gtf)
gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]

# Exons: GRanges named by transcript_id, with exon_position metadata
st_exons_raw  <- gtf_clean[gtf_clean$type == "exon"]

# Build position-aware exon table first (needed for intron derivation too)
# Terminal exon position (first/last/single) is flagged here and carried into
# st_exons_flat as metadata for use in the section 5a exonic filter:
#   - "internal" exons: must be fully contained within a reference exon block
#   - "first" exons (lowest genomic coord): only their END (splice-site side)
#     must land within a reference exon — the 5' terminus is free to extend
#   - "last" exons (highest genomic coord): only their START (splice-site side)
#     must land within a reference exon — the 3' terminus is free to extend
#   - "single" exon transcripts have no splice sites → skip the exonic filter
# Position is defined in genomic (not transcript) coordinates and is therefore
# strand-agnostic: first = lowest start, last = highest end.
exon_dt <- data.table(
    chr           = as.character(seqnames(st_exons_raw)),
    start         = start(st_exons_raw),
    end           = end(st_exons_raw),
    strand        = as.character(strand(st_exons_raw)),
    transcript_id = st_exons_raw$transcript_id
)[order(transcript_id, start)]

exon_dt[, exon_position := {
    n <- .N
    if (n == 1L) "single"
    else c("first", rep("internal", n - 2L), "last")
}, by = transcript_id]

intron_dt <- exon_dt[, {
    n <- .N
    if (n < 2L) {
        data.table(chr = character(), intron_start = integer(),
                   intron_end = integer(), strand = character())
    } else {
        data.table(
            chr          = chr[seq_len(n - 1L)],
            intron_start = end[seq_len(n - 1L)] + 1L,
            intron_end   = start[2L:n] - 1L,
            strand       = strand[seq_len(n - 1L)]
        )
    }
}, by = transcript_id]
intron_dt <- intron_dt[intron_start <= intron_end]

# Build st_exons_flat now that exon_dt (with exon_position) is available
st_exons_flat <- st_exons_raw
names(st_exons_flat) <- st_exons_flat$transcript_id
# Match exon_position from exon_dt by transcript_id + start coordinate
ep_lookup <- exon_dt[, .(transcript_id, start, exon_position)]
setkey(ep_lookup, transcript_id, start)
st_exons_flat$exon_position <- ep_lookup[
    J(st_exons_flat$transcript_id, start(st_exons_flat)),
    exon_position
]

st_introns_flat <- GRanges(
    seqnames = intron_dt$chr,
    ranges   = IRanges(start = intron_dt$intron_start, end = intron_dt$intron_end),
    strand   = intron_dt$strand
)
names(st_introns_flat) <- intron_dt$transcript_id

# Rebuild GRL equivalents (needed for sections 9 and for novel tx detection)
st_exons_grl   <- split(st_exons_flat,   names(st_exons_flat))
st_introns_grl <- split(st_introns_flat, names(st_introns_flat))

# Novel transcripts: MSTRG prefix from StringTie denotes non-reference assemblies
novel_st_tx_names <- unique(c(
    names(st_introns_grl)[grepl("^MSTRG", names(st_introns_grl))],
    names(st_exons_grl)[grepl("^MSTRG",  names(st_exons_grl))]
))
cat(sprintf("STRINGTIE: %d novel (MSTRG) transcripts identified.\n", length(novel_st_tx_names)))

# ── 5. NEGATIVE FILTER: REMOVE TRANSCRIPTS WITH UNSUPPORTED NOVEL STRUCTURE ───
# Two classes of structural novelty are evaluated per transcript:
#
# (a) NOVEL EXONIC CONTIGS: portions of MSTRG exons not covered by any reference
#     exon. Each contig >= min_aad_length bp must be overlapped by at least one
#     positive dPSI Whippet event (CE+, AA+, AD+, RI+).
#
# (b) NOVEL INTRONIC CONTIGS: portions of MSTRG introns that overlap reference
#     exons (reference-exonic sequence spliced out in the novel transcript). Each
#     such contig >= min_aad_length bp must be overlapped by at least one negative
#     dPSI Whippet event (CE-, AA-, AD-, RI-).
#
# IMPLEMENTATION: vectorised over the flat GRanges with tx_id carried as names,
# avoiding both the genome-wide reduce() pooling problem (which loses per-
# transcript identity) and the per-transcript Filter() loop overhead.
#
# Section 5 uses SIMPLE OVERLAP for Whippet confirmation — boundary-aware logic
# belongs only in the positive filter (section 6).

cat("FILTER: Applying negative filter (unsupported novel exonic/intronic contigs)...\n")

# All positive and negative dPSI Whippet events as flat GRanges
whippet_pos <- whippet_gr[whippet_gr$dpsi > 0]
whippet_neg <- whippet_gr[whippet_gr$dpsi < 0]

# Restrict MSTRG introns only — names carry tx_id throughout
novel_st_introns_flat <- st_introns_flat[names(st_introns_flat) %in% novel_st_tx_names]

# ── 5a. Novel exonic contigs ──────────────────────────────────────────────────
# Exons are handled differently by position:
#
#   SINGLE-EXON transcripts: no splice sites exist on either end, so there is
#     no constrained boundary. These are skipped entirely.
#
#   INTERNAL exons: both boundaries are splice sites. The full exon must be
#     contained within a reference exon block — any extension is novel content.
#
#   FIRST exons (lowest genomic coord): the END is the splice-site boundary and
#     must land within a reference exon (within boundary_tol). The START is the
#     transcript terminus and is free to extend — UTR extension is expected noise.
#
#   LAST exons (highest genomic coord): the START is the splice-site boundary and
#     must land within a reference exon (within boundary_tol). The END is free.
#
# In all cases, single-nucleotide wobble up to boundary_tol is permitted on the
# constrained boundary before flagging.

novel_st_exons_flat <- st_exons_flat[names(st_exons_flat) %in% novel_st_tx_names]

exons_single   <- novel_st_exons_flat[novel_st_exons_flat$exon_position == "single"]
exons_internal <- novel_st_exons_flat[novel_st_exons_flat$exon_position == "internal"]
exons_first    <- novel_st_exons_flat[novel_st_exons_flat$exon_position == "first"]
exons_last     <- novel_st_exons_flat[novel_st_exons_flat$exon_position == "last"]

# Internal exons: both START and END are splice-site boundaries.
# Check each boundary independently — it must land within a reference exon block
# within boundary_tol. This tolerates minor splice-site coordinate wobble from
# StringTie assembly (typically 1–3bp) without falsely flagging the transcript.
# The previous type="within" check required the ENTIRE exon to fit inside one
# reference block, which failed for any exon with boundary wobble even when
# both splice sites are canonical.
check_internal <- function(exons) {
    if (length(exons) == 0L) return(character(0))
    # 1-bp GRanges at start boundary and at end boundary of each exon
    bnd_start <- GRanges(seqnames = seqnames(exons),
                         ranges   = IRanges(start = start(exons), width = 1L),
                         strand   = strand(exons))
    bnd_end   <- GRanges(seqnames = seqnames(exons),
                         ranges   = IRanges(start = end(exons), width = 1L),
                         strand   = strand(exons))
    # Both boundaries must land within a ref exon block (within boundary_tol)
    ok_start <- unique(queryHits(findOverlaps(bnd_start, ref_exons_reduced,
                                              maxgap = boundary_tol,
                                              ignore.strand = FALSE)))
    ok_end   <- unique(queryHits(findOverlaps(bnd_end, ref_exons_reduced,
                                              maxgap = boundary_tol,
                                              ignore.strand = FALSE)))
    both_ok  <- intersect(ok_start, ok_end)
    not_ok   <- setdiff(seq_along(exons), both_ok)
    if (length(not_ok) == 0L) return(character(0))
    # Remaining: check Whippet+ confirmation (simple overlap with full exon range)
    unconf <- exons[not_ok]
    if (length(whippet_pos) > 0L) {
        conf <- unique(queryHits(findOverlaps(unconf, whippet_pos,
                                             ignore.strand = FALSE)))
        if (length(conf) > 0L) unconf <- unconf[-conf]
    }
    # Only flag if the unconfirmed exon is substantially novel (>= min_aad_length
    # extension beyond any reference block on at least one side)
    if (length(unconf) == 0L) return(character(0))
    any_ref <- findOverlaps(unconf, ref_exons_reduced, ignore.strand = FALSE)
    no_ref  <- setdiff(seq_along(unconf), unique(queryHits(any_ref)))
    sig <- names(unconf[no_ref[width(unconf[no_ref]) >= min_aad_length]])
    if (length(unique(queryHits(any_ref))) > 0L) {
        q   <- unconf[queryHits(any_ref)]
        s   <- ref_exons_reduced[subjectHits(any_ref)]
        ext <- tapply(pmax(pmax(0L, start(s) - start(q)),
                           pmax(0L, end(q)   - end(s))),
                      queryHits(any_ref), max)
        sig <- unique(c(sig,
                        names(unconf)[as.integer(names(ext))[ext >= min_aad_length]]))
    }
    sig
}
sig_internal <- check_internal(exons_internal)

# Helper: check that a single boundary of a terminal exon lands within a ref exon.
# `coord_fn` extracts the constrained coordinate (start or end) from each exon.
check_terminal <- function(exons, coord_fn) {
    if (length(exons) == 0L) return(character(0))
    # Build a 1-bp GRanges at the splice-site boundary of each exon
    bnd <- GRanges(seqnames = seqnames(exons),
                   ranges   = IRanges(start = coord_fn(exons),
                                      width = 1L),
                   strand   = strand(exons))
    names(bnd) <- names(exons)
    # Boundary lands in a ref exon block (within tolerance) → no novel content
    hits_ref <- findOverlaps(bnd, ref_exons_reduced,
                             maxgap = boundary_tol, ignore.strand = FALSE)
    ref_ok   <- unique(queryHits(hits_ref))
    # Boundaries not in any ref exon: check Whippet+ confirmation
    not_ok <- setdiff(seq_along(bnd), ref_ok)
    if (length(not_ok) == 0L) return(character(0))
    unconf <- bnd[not_ok]
    if (length(whippet_pos) > 0L) {
        conf <- unique(queryHits(findOverlaps(exons[not_ok], whippet_pos,
                                             ignore.strand = FALSE)))
        unconf <- unconf[-conf]
    }
    unique(names(unconf))
}

sig_first <- check_terminal(exons_first, end)    # END is splice-site side
sig_last  <- check_terminal(exons_last,  start)  # START is splice-site side

tx_unsupported_exonic <- unique(c(sig_internal, sig_first, sig_last))
cat(sprintf("FILTER: %d transcripts carry unsupported novel exonic contigs.\n",
            length(tx_unsupported_exonic)))

# ── 5b. Novel intronic contigs ────────────────────────────────────────────────
# CRITICAL: pre-clear MSTRG introns that exactly match a reference intron.
# Any intron with a reference counterpart represents canonical splicing — it
# correctly splices out reference-exonic territory from an alternative isoform
# and needs no Whippet confirmation. Skipping this step causes reduce()-merged
# reference exon coverage to falsely flag constitutive introns as novel, which
# in multi-isoform genes eliminates nearly all transcripts.
#
# Only truly novel introns — those with no reference match within boundary_tol —
# are then checked against reference exon coverage and Whippet- events.
if (length(novel_st_introns_flat) > 0L) {
    ref_match     <- findOverlaps(novel_st_introns_flat, ref_introns_gr,
                                  type = "equal", maxgap = boundary_tol,
                                  ignore.strand = FALSE)
    canonical_idx <- unique(queryHits(ref_match))
    truly_novel_introns <- if (length(canonical_idx) > 0L)
        novel_st_introns_flat[-canonical_idx] else novel_st_introns_flat

    cat(sprintf("FILTER: %d / %d MSTRG introns cleared as canonical (reference match).\n",
                length(canonical_idx), length(novel_st_introns_flat)))

    if (length(truly_novel_introns) > 0L) {
        ov_int <- findOverlaps(truly_novel_introns, ref_exons_reduced, ignore.strand = FALSE)

        if (length(ov_int) > 0L) {
            intron_sub <- truly_novel_introns[queryHits(ov_int)]
            ref_sub    <- ref_exons_reduced[subjectHits(ov_int)]
            pint_all   <- pintersect(intron_sub, ref_sub, ignore.strand = TRUE)
            keep_idx   <- which(width(pint_all) >= min_aad_length)
            novel_intronic        <- pint_all[keep_idx]
            names(novel_intronic) <- names(intron_sub)[keep_idx]

            if (length(novel_intronic) > 0L && length(whippet_neg) > 0L) {
                conf_in_hits  <- findOverlaps(novel_intronic, whippet_neg, ignore.strand = FALSE)
                confirmed_in  <- unique(queryHits(conf_in_hits))
                unconf_intron <- if (length(confirmed_in) > 0L)
                    novel_intronic[-confirmed_in] else novel_intronic
            } else {
                unconf_intron <- novel_intronic
            }
            tx_unsupported_intronic <- unique(names(unconf_intron))
        } else {
            tx_unsupported_intronic <- character(0)
        }
    } else {
        tx_unsupported_intronic <- character(0)
    }
} else {
    tx_unsupported_intronic <- character(0)
}
cat(sprintf("FILTER: %d transcripts carry unsupported novel intronic contigs.\n",
            length(tx_unsupported_intronic)))

unsupported_tx <- unique(c(tx_unsupported_exonic, tx_unsupported_intronic))
cat(sprintf("FILTER: %d transcripts flagged in total.\n", length(unsupported_tx)))

# ── 6. POSITIVE FILTER: HELPERS ───────────────────────────────────────────────
# Each helper returns the names of StringTie transcripts that support a given
# class of Whippet event. The boundary_tol applies to the "containment" side;
# for AA/AD the CE-adjacent boundary uses the same tolerance on the stricter
# side (see inline comments for strand-aware logic).

# CE/RI positive dPSI — StringTie exon fully contains the node
tx_exon_contains <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    ok   <- start(ex_h) <= start(ev_h) + boundary_tol &
            end(ex_h)   >= end(ev_h)   - boundary_tol
    unique(names(ex_h)[ok])
}

# CE/RI negative dPSI — StringTie intron fully spans the node
tx_intron_spans <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    ok   <- start(in_h) <= start(ev_h) + boundary_tol &
            end(in_h)   >= end(ev_h)   - boundary_tol
    unique(names(in_h)[ok])
}

# AA positive dPSI — exon spans node; alternative acceptor boundary matches.
# The alternative acceptor is the "far" boundary from CE:
#   + strand → LEFT boundary of node  (lower genomic coord)
#   - strand → RIGHT boundary of node (higher genomic coord)
# A downstream exon cannot satisfy this: it begins AFTER the node ends (+),
# so it cannot span the node at all.
tx_aa_pos <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    # boundary_tol applied to both sides of spans: StringTie wobble at the
    # CE-adjacent (right on +) boundary is as likely as at the alt-acceptor side
    spans <- start(ex_h) <= start(ev_h) + boundary_tol &
             end(ex_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+",
                    abs(start(ex_h) - start(ev_h)) <= boundary_tol,  # alt acceptor = left on +
                    abs(end(ex_h)   - end(ev_h))   <= boundary_tol)  # alt acceptor = right on -
    unique(names(ex_h)[spans & bnd])
}

# AA negative dPSI — intron spans node; canonical acceptor boundary matches.
# The canonical acceptor is the CE-adjacent boundary:
#   + strand → RIGHT boundary of node
#   - strand → LEFT boundary of node
tx_aa_neg <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    spans <- start(in_h) <= start(ev_h) + boundary_tol &
             end(in_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+",
                    abs(end(in_h)   - end(ev_h))   <= boundary_tol,  # canonical acceptor = right on +
                    abs(start(in_h) - start(ev_h)) <= boundary_tol)  # canonical acceptor = left on -
    unique(names(in_h)[spans & bnd])
}

# AD positive dPSI — exon spans node; alternative donor boundary matches.
# The alternative donor is the "far" boundary from CE:
#   + strand → RIGHT boundary of node
#   - strand → LEFT boundary of node
# An upstream exon cannot satisfy this: it ends BEFORE the node starts (+),
# so it cannot span the node at all.
tx_ad_pos <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    spans <- start(ex_h) <= start(ev_h) + boundary_tol &
             end(ex_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+",
                    abs(end(ex_h)   - end(ev_h))   <= boundary_tol,  # alt donor = right on +
                    abs(start(ex_h) - start(ev_h)) <= boundary_tol)  # alt donor = left on -
    unique(names(ex_h)[spans & bnd])
}

# AD negative dPSI — intron spans node; canonical donor boundary matches.
# The canonical donor is the CE-adjacent boundary:
#   + strand → LEFT boundary of node
#   - strand → RIGHT boundary of node
tx_ad_neg <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]
    ev_h <- events[subjectHits(hits)]
    spans <- start(in_h) <= start(ev_h) + boundary_tol &
             end(in_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+",
                    abs(start(in_h) - start(ev_h)) <= boundary_tol,  # canonical donor = left on +
                    abs(end(in_h)   - end(ev_h))   <= boundary_tol)  # canonical donor = right on -
    unique(names(in_h)[spans & bnd])
}

# ── 7. POSITIVE FILTER: APPLY ────────────────────────────────────────────────
cat("FILTER: Applying positive filter (StringTie support for Whippet events)...\n")

supported_by_whippet <- unique(c(
    tx_exon_contains(ce_pos, st_exons_flat),
    tx_intron_spans( ce_neg, st_introns_flat),
    tx_exon_contains(ri_pos, st_exons_flat),
    tx_intron_spans( ri_neg, st_introns_flat),
    tx_aa_pos(aa_pos, st_exons_flat),
    tx_aa_neg(aa_neg, st_introns_flat),
    tx_ad_pos(ad_pos, st_exons_flat),
    tx_ad_neg(ad_neg, st_introns_flat)
))

# ── 8. COMBINE FILTERS ────────────────────────────────────────────────────────
valid_tx_names <- intersect(
    supported_by_whippet,
    setdiff(novel_st_tx_names, unsupported_tx)
)
cat(sprintf("FILTER: %d novel transcripts pass both filters.\n", length(valid_tx_names)))

# ── 9. IDENTIFY CONFIRMED WHIPPET EVENTS ─────────────────────────────────────
# Re-run the same overlap logic restricted to valid transcripts only, returning
# lsv_ids rather than transcript names, to identify which events are confirmed.
cat("EXPORT: Identifying confirmed Whippet events...\n")

st_final_introns    <- st_introns_flat[names(st_introns_flat) %in% valid_tx_names]
st_final_exons_flat <- unlist(st_exons_grl[names(st_exons_grl) %in% valid_tx_names])

# Parallel helpers — same logic as section 6 but return lsv_ids of confirmed events
confirm_exon_contains <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    ok   <- start(ex_h) <= start(ev_h) + boundary_tol & end(ex_h) >= end(ev_h) - boundary_tol
    unique(ev_h$lsv_id[ok])
}

confirm_intron_spans <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    ok   <- start(in_h) <= start(ev_h) + boundary_tol & end(in_h) >= end(ev_h) - boundary_tol
    unique(ev_h$lsv_id[ok])
}

confirm_aa_pos <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    spans <- start(ex_h) <= start(ev_h) + boundary_tol &
             end(ex_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+", abs(start(ex_h) - start(ev_h)) <= boundary_tol,
                                  abs(end(ex_h)   - end(ev_h))   <= boundary_tol)
    unique(ev_h$lsv_id[spans & bnd])
}

confirm_aa_neg <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    spans <- start(in_h) <= start(ev_h) + boundary_tol &
             end(in_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+", abs(end(in_h)   - end(ev_h))   <= boundary_tol,
                                  abs(start(in_h) - start(ev_h)) <= boundary_tol)
    unique(ev_h$lsv_id[spans & bnd])
}

confirm_ad_pos <- function(events, exons_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(exons_flat, events, type = "any", ignore.strand = FALSE)
    ex_h <- exons_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    spans <- start(ex_h) <= start(ev_h) + boundary_tol &
             end(ex_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+", abs(end(ex_h)   - end(ev_h))   <= boundary_tol,
                                  abs(start(ex_h) - start(ev_h)) <= boundary_tol)
    unique(ev_h$lsv_id[spans & bnd])
}

confirm_ad_neg <- function(events, introns_flat) {
    if (length(events) == 0L) return(character(0))
    events <- events[width(events) >= min_aad_length]
    if (length(events) == 0L) return(character(0))
    hits <- findOverlaps(introns_flat, events, type = "any", ignore.strand = FALSE)
    in_h <- introns_flat[queryHits(hits)]; ev_h <- events[subjectHits(hits)]
    spans <- start(in_h) <= start(ev_h) + boundary_tol &
             end(in_h)   >= end(ev_h)   - boundary_tol
    str_h <- as.character(strand(ev_h))
    bnd   <- ifelse(str_h == "+", abs(start(in_h) - start(ev_h)) <= boundary_tol,
                                  abs(end(in_h)   - end(ev_h))   <= boundary_tol)
    unique(ev_h$lsv_id[spans & bnd])
}

confirmed_lsv_ids <- unique(c(
    confirm_exon_contains(ce_pos, st_final_exons_flat),
    confirm_intron_spans( ce_neg, st_final_introns),
    confirm_exon_contains(ri_pos, st_final_exons_flat),
    confirm_intron_spans( ri_neg, st_final_introns),
    confirm_aa_pos(aa_pos, st_final_exons_flat),
    confirm_aa_neg(aa_neg, st_final_introns),
    confirm_ad_pos(ad_pos, st_final_exons_flat),
    confirm_ad_neg(ad_neg, st_final_introns)
))

# ── 10. EXPORT ────────────────────────────────────────────────────────────────
final_whippet_events <- whippet_events %>% filter(lsv_id %in% confirmed_lsv_ids)

cat(sprintf("EXPORT: %d / %d Whippet events confirmed by novel StringTie transcripts.\n",
            nrow(final_whippet_events), nrow(whippet_events)))

outfile <- file.path(output_dir, "whippet_confirmed_events.tsv")
fwrite(final_whippet_events, outfile, sep = "\t")
cat(sprintf("EXPORT: Confirmed events written to %s\n", outfile))

cat(sprintf("EXPORT: Generating validated consensus matrix (supported by StringTie)..."))
gtf_novel_only <- gtf_clean[mcols(gtf_clean)$transcript_id %in% valid_tx_names]
rtracklayer::export(gtf_novel_only, novel_gtf, format="gtf")
cat(sprintf("EXPORT: Novel-only GTF Construction Complete. File Saved at %s",novel_gtf))