import sys
from pathlib import Path

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

REQUIRED_COLS = {"chr", "js", "je", "strand"}


def load_df(path: str, label: str = None, tool_name: str = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    # If necessary, normalize start/end to js/je
    if "start" in df.columns and "end" in df.columns:
        if "js" not in df.columns:
            df.rename(columns={"start": "js"}, inplace=True)
        if "je" not in df.columns:
            df.rename(columns={"end": "je"}, inplace=True)

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        label_str = f" ({label})" if label else ""
        sys.exit(f"[ERROR]{label_str} ({path}) is missing columns: {missing}")
    df["js"] = df["js"].astype(int)
    df["je"] = df["je"].astype(int)
    df["strand"] = df["strand"].astype(str)
    df["chr"] = df["chr"].astype(str)
    return df


# ---------------------------------------------------------------------------
# Matching logic
# ---------------------------------------------------------------------------

def match_junctions(
    truth: pd.DataFrame,
    pred: pd.DataFrame,
    tol: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match every truth junction to the best available predicted junction.

    Priority: exact match > partial match > fn (no match).
    Each predicted junction is used as the *primary* match for at most one
    truth junction (greedy, truth-row order). Unmatched predictions are FPs.
    """
    # Build per-chr+strand lookup for speed
    pred_by_chr_strand: dict = {}
    for idx, row in pred.iterrows():
        key = (row["chr"], row["strand"])
        pred_by_chr_strand.setdefault(key, []).append(
            (idx, int(row["js"]), int(row["je"]))
        )

    truth_rows   = []
    pred_matched = set()

    for _, trow in truth.iterrows():
        tjs, tje = int(trow["js"]), int(trow["je"])

        key        = (trow["chr"], trow["strand"])
        candidates = pred_by_chr_strand.get(key, [])

        best_status = "fn"
        best_pjs, best_pje = np.nan, np.nan
        best_pidx   = None

        for pidx, pjs, pje in candidates:
            js_match = abs(tjs - pjs) <= tol
            je_match = abs(tje - pje) <= tol

            if js_match and je_match:
                best_status = "exact"
                best_pjs, best_pje = pjs, pje
                best_pidx = pidx
                break                        # exact is the best possible
            elif js_match or je_match:
                if best_status != "exact":   # don't downgrade from exact
                    best_status = "partial"
                    best_pjs, best_pje = pjs, pje
                    best_pidx = pidx

        if best_pidx is not None:
            pred_matched.add(best_pidx)

        truth_rows.append({
            **trow.to_dict(),
            "match_status":    best_status,
            "matched_pred_js": best_pjs,
            "matched_pred_je": best_pje,
        })

    # Label predictions
    pred_res = pred.copy()
    pred_res["match_status"] = pred_res.index.map(
        lambda i: "matched" if i in pred_matched else "fp"
    )

    return pd.DataFrame(truth_rows), pred_res


def match_contigs_exon_based(
    truth: pd.DataFrame,
    pred: pd.DataFrame,
    tol: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Exon contig matching with boundary matching only.

    For each truth contig, find best matching predicted contig:
      exact   - both boundaries match within tolerance
      partial - only one boundary matches within tolerance
      fn      - neither boundary matches

    Each predicted contig used as primary match for at most one truth contig.
    """
    pred_by_chr_strand: dict = {}
    for idx, row in pred.iterrows():
        key = (row["chr"], row["strand"])
        pred_by_chr_strand.setdefault(key, []).append(
            (idx, int(row["js"]), int(row["je"]))
        )

    truth_rows   = []
    pred_matched = set()

    for _, trow in truth.iterrows():
        tjs, tje = int(trow["js"]), int(trow["je"])

        key        = (trow["chr"], trow["strand"])
        candidates = pred_by_chr_strand.get(key, [])

        best_status = "fn"
        best_pjs, best_pje = np.nan, np.nan
        best_pidx   = None

        for pidx, pjs, pje in candidates:
            js_match = abs(tjs - pjs) <= tol
            je_match = abs(tje - pje) <= tol

            if js_match and je_match:
                best_status = "exact"
                best_pjs, best_pje = pjs, pje
                best_pidx = pidx
                break
            elif js_match or je_match:
                if best_status != "exact":
                    best_status = "partial"
                    best_pjs, best_pje = pjs, pje
                    best_pidx = pidx

        if best_pidx is not None:
            pred_matched.add(best_pidx)

        truth_rows.append({
            **trow.to_dict(),
            "match_status":    best_status,
            "matched_pred_js": best_pjs,
            "matched_pred_je": best_pje,
        })

    pred_res = pred.copy()
    pred_res["match_status"] = pred_res.index.map(
        lambda i: "matched" if i in pred_matched else "fp"
    )

    return pd.DataFrame(truth_rows), pred_res


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def _r4(x) -> float:
    """Round to 4 decimal places; return NaN if not finite."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    return round(float(x), 4)


def compute_metrics(truth_res: pd.DataFrame, pred_res: pd.DataFrame) -> dict:
    """
    Compute boundary-based scalar metrics for one (tool, event_type) slice.

    truth_res  - annotated truth junctions (possibly filtered to one event_type)
    pred_res   - annotated predictions for the whole tool (FP count is
                 tool-wide; event_type is not available for predictions)
    """
    n_truth = len(truth_res)
    exact   = int((truth_res["match_status"] == "exact").sum())
    partial = int((truth_res["match_status"] == "partial").sum())
    fn      = int((truth_res["match_status"] == "fn").sum())
    fp      = int((pred_res["match_status"]  == "fp").sum())
    n_pred  = len(pred_res)

    detected     = exact + partial
    precision    = detected / n_pred  if n_pred  > 0 else np.nan
    recall_exact = exact    / n_truth if n_truth > 0 else np.nan
    recall_any   = detected / n_truth if n_truth > 0 else np.nan

    def f1(p, r):
        return 2 * p * r / (p + r) if (not np.isnan(p)) and (p + r) > 0 else np.nan

    return {
        "n_truth":              n_truth,
        "n_pred":               n_pred,
        "exact":                exact,
        "partial":              partial,
        "fn":                   fn,
        "fp":                   fp,
        "detection_rate_exact": _r4(recall_exact),
        "detection_rate_any":   _r4(recall_any),
        "precision":            _r4(precision),
        "f1_exact":             _r4(f1(precision, recall_exact)),
        "f1_any":               _r4(f1(precision, recall_any)),
    }


def compute_metrics_contigs_exon_based(truth_res: pd.DataFrame, pred_res: pd.DataFrame) -> dict:
    """Compute boundary-based metrics for exon contigs."""
    n_truth = len(truth_res)
    exact   = int((truth_res["match_status"] == "exact").sum())
    partial = int((truth_res["match_status"] == "partial").sum())
    fn      = int((truth_res["match_status"] == "fn").sum())
    fp      = int((pred_res["match_status"]  == "fp").sum())
    n_pred  = len(pred_res)

    detected     = exact + partial
    precision    = detected / n_pred  if n_pred  > 0 else np.nan
    recall_exact = exact    / n_truth if n_truth > 0 else np.nan
    recall_any   = detected / n_truth if n_truth > 0 else np.nan

    def f1(p, r):
        return 2 * p * r / (p + r) if (not np.isnan(p)) and (p + r) > 0 else np.nan

    return {
        "n_truth":              n_truth,
        "n_pred":               n_pred,
        "exact":                exact,
        "partial":              partial,
        "fn":                   fn,
        "fp":                   fp,
        "detection_rate_exact": _r4(recall_exact),
        "detection_rate_any":   _r4(recall_any),
        "precision":            _r4(precision),
        "f1_exact":             _r4(f1(precision, recall_exact)),
        "f1_any":               _r4(f1(precision, recall_any)),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Load both ground-truth files via named inputs (not positional indices).
    # Exon-based tools are scored against exon contigs; all other tools against junctions.
    truth_junctions = load_df(snakemake.input.ground_truth_junctions, "ground truth junctions")
    truth_exons     = load_df(snakemake.input.ground_truth_exons,     "ground truth exons",
                              tool_name="exon_tool")  # normalises start/end -> js/je

    # Assign a shared integer event_id (row order) to both truth files.
    # Junction and exon truth files describe the same biological events in the
    # same order but use different coordinate systems (intron vs. exon bounds).
    # Using event_id as the pivot key ensures that junction tools and exon tools
    # contribute to the same row in the overview tables.
    truth_junctions = truth_junctions.sort_values(by=["chr", "js"]).reset_index(drop=True)
    truth_exons     = truth_exons.sort_values(by=["chr", "js"]).reset_index(drop=True)
    truth_junctions["event_id"] = truth_junctions.index
    truth_exons["event_id"]     = truth_exons.index

    out_dir = Path(snakemake.output[0]).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows        = []
    per_junction_frames = []
    fp_frames           = []

    pred_paths = list(snakemake.input.tool_outputs)

    # exon_tools: set of tool names that output exon contigs rather than splice
    # junctions (e.g. Whippet, GTF_Whippet). Populated from snakemake.params so
    # that adding a new exon-based tool only requires updating the Snakefile /
    # config — no changes to this script are needed.
    exon_tools = set(snakemake.params.exon_tools)

    for tool_name, pred_path in zip(snakemake.params.tool_names, pred_paths):
        is_exon_tool = tool_name in exon_tools

        # Select the appropriate ground truth for this tool
        truth          = truth_exons     if is_exon_tool else truth_junctions
        has_event_type = "event_type" in truth.columns
        event_types    = truth["event_type"].unique() if has_event_type else []

        pred = load_df(pred_path, tool_name, tool_name=tool_name)

        if is_exon_tool:
            truth_res, pred_res = match_contigs_exon_based(truth, pred, snakemake.params.tolerance)
            m = compute_metrics_contigs_exon_based(truth_res, pred_res)
        else:
            truth_res, pred_res = match_junctions(truth, pred, snakemake.params.tolerance)
            m = compute_metrics(truth_res, pred_res)

        summary_rows.append({"tool": tool_name, "event_type": "ALL", **m})

        # per-event_type metrics
        if has_event_type:
            for etype in event_types:
                t_sub = truth_res[truth_res["event_type"] == etype]
                if is_exon_tool:
                    m_et = compute_metrics_contigs_exon_based(t_sub, pred_res)
                else:
                    m_et = compute_metrics(t_sub, pred_res)
                summary_rows.append({"tool": tool_name, "event_type": str(etype), **m_et})

        truth_res["tool"] = tool_name
        per_junction_frames.append(truth_res)

        fp_df = pred_res[pred_res["match_status"] == "fp"].copy()
        fp_df["tool"] = tool_name
        fp_frames.append(fp_df)

    # ---- assemble and write outputs ----

    summary_df = pd.DataFrame(summary_rows)
    pj_df      = pd.concat(per_junction_frames, ignore_index=True)
    fp_df_all  = pd.concat(fp_frames, ignore_index=True) if fp_frames else pd.DataFrame()

    ## summary df
    col_order = [
        "tool", "event_type",
        "n_truth", "n_pred",
        "exact", "partial", "fn", "fp",
        "detection_rate_exact", "detection_rate_any",
        "precision", "f1_exact", "f1_any",
    ]
    summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]

    # ---- overview table: 1 row per event_type, 1 column per tool ----
    # Each cell: "exact=X | partial=Y | fn=Z"
    overview_rows = summary_df[summary_df["event_type"] != "ALL"].copy()
    overview_rows["cell"] = (
        "exact="   + overview_rows["exact"].astype(str) + " | " +
        "partial=" + overview_rows["partial"].astype(str) + " | " +
        "fn="      + overview_rows["fn"].astype(str)
    )
    overview_df = overview_rows.pivot(index="event_type", columns="tool", values="cell")
    overview_df.columns.name = None
    overview_df.index.name   = "event_type"
    # Guarantee consistent tool column order
    tool_order = [t for t in snakemake.params.tool_names if t in overview_df.columns]
    overview_df = overview_df[tool_order]

    # ---- per-junction overview: 1 row per splice event, 1 column per tool ----
    # Cell content: match_status.
    pj_df["cell"] = pj_df["match_status"]
    junction_key = ["event_id", "event_type"]
    pj_overview_df = pj_df.pivot(index=junction_key, columns="tool", values="cell")
    pj_overview_df.columns.name = None
    pj_overview_df = pj_overview_df.reset_index()
    # Consistent tool column order
    meta_cols = junction_key
    tool_cols  = [t for t in snakemake.params.tool_names if t in pj_overview_df.columns]
    pj_overview_df = pj_overview_df[meta_cols + tool_cols]
    # Sort by event_type then event_id for readability
    pj_overview_df = pj_overview_df.sort_values(["event_type", "event_id"])

    # ---- write outputs ----
    summary_path = snakemake.output[0]
    out_dir = Path(summary_path).parent
    prefix = Path(summary_path).stem.rsplit("_summary", 1)[0]

    pj_path       = str(out_dir / f"{prefix}_per_junction.tsv")
    fp_path       = str(out_dir / f"{prefix}_fp.tsv")
    overview_path = str(out_dir / f"{prefix}_overview.tsv")
    pj_overview_path = str(out_dir / f"{prefix}_junction_overview.tsv")

    summary_df.to_csv(summary_path,      sep="\t", index=False)
    pj_df.to_csv(pj_path,                sep="\t", index=False)
    fp_df_all.to_csv(fp_path,            sep="\t", index=False)
    overview_df.to_csv(overview_path,    sep="\t")
    pj_overview_df.to_csv(pj_overview_path, sep="\t", index=False)

    # ---- console summary ----
    print(f"\nTolerance : +/-{snakemake.params.tolerance} bp")
    print(f"Tools     : {', '.join(snakemake.params.tool_names)}\n")

    display_cols = [
        "tool", "event_type",
        "n_truth", "exact", "partial", "fn", "fp",
        "detection_rate_any", "precision",
    ]
    display_cols = [c for c in display_cols if c in summary_df.columns]
    all_summary  = summary_df[summary_df["event_type"] == "ALL"]
    print(all_summary[display_cols].to_string(index=False))

    print(f"\nOutputs written:")
    print(f"  {summary_path:<55}  <- per-tool aggregate metrics")
    print(f"  {pj_path:<55}  <- every truth junction annotated")
    print(f"  {fp_path:<55}  <- false positives per tool")
    print(f"  {overview_path:<55}  <- event_type x tool overview")
    print(f"  {pj_overview_path:<55}  <- per-junction x tool overview")


if __name__ == "__main__":
    main()