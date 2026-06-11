## Preparation of tool output for scoring against ground truth. This script is borderline useless and 
## could without much effort be integrated into the following scoring script. (The only "preparation" that
## happening is renaming the start/end columns to js/je).
## Though since it does no harm and conveniently renames all output files to the same naming 
## convention used by all the other snakmake scripts, I have kept it so far.
## During the time it took me to write this I could have probably done half of the integration already.
## If you're still reading, consider yourself ragebaited.

cat("
Prepare output to be scored
Usage: Rscript prepare_output.R <tool_output_file> <prepared_output_file>
")
suppressPackageStartupMessages({
    library(dplyr) 
    library(stringr)
    library(readr)
    library(tidyr)
})
log_file <- snakemake@log[[1]]
# expected input: tool output with columns chr, start, end, strand
# Expect a named parameter `cfg` in snakemake@params (see Snakefile)
tool_outputs <- NULL
if (!is.null(snakemake@params$cfg)) {
    tool_outputs <- snakemake@params$cfg
} else if (length(snakemake@params) >= 1) {
    tool_outputs <- snakemake@params[[1]]
}

# Outputs (one file per tool)
output_files <- snakemake@output

if (is.null(tool_outputs) || length(tool_outputs) == 0) {
    stop("No tool output files provided in snakemake params (cfg)")
}

if (length(tool_outputs) != length(output_files)) {
    warning(sprintf("Number of tool inputs (%d) != number of outputs (%d). Proceeding with min length.",
                    length(tool_outputs), length(output_files)))
}

N <- min(length(tool_outputs), length(output_files))

for (i in seq_len(N)) {
    tool_output_file <- tool_outputs[[i]]
    output_file <- output_files[[i]]

    cat(sprintf("Processing tool output: %s -> %s\n", tool_output_file, output_file))

    # Read input, be permissive about missing cols
    tool_output <- tryCatch(
        read.table(tool_output_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, na.strings = c("", "NA")),
        error = function(e) stop(sprintf("Failed to read tool output '%s': %s", tool_output_file, conditionMessage(e)))
    )

    # Ensure required columns exist or can be renamed
    if (!all(c("chr", "js", "je", "strand") %in% colnames(tool_output))) {
        if (all(c("chr", "start", "end", "strand") %in% colnames(tool_output))) {
            cat("Renaming columns 'start'/'end' -> 'js'/'je'\n")
            colnames(tool_output)[colnames(tool_output) == "start"] <- "js"
            colnames(tool_output)[colnames(tool_output) == "end"] <- "je"
        } else {
            stop(sprintf("Input file '%s' must contain columns: chr,start,end,strand OR chr,js,je,strand", tool_output_file))
        }
    }

    # Select the columns we need
    prepared_output <- tool_output %>%
        dplyr::select(chr, js, je, strand)

    # Ensure parent directory exists
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir) && out_dir != ".") {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        cat(sprintf("Created directory %s\n", out_dir))
    }

    # Write prepared output
    write.table(prepared_output, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    cat(sprintf("Wrote prepared output to %s\n", output_file))
}
