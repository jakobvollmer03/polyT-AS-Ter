#!/bin/bash

## Prepare transcript sequences from StringTie GTF for Polyester simulation
## This script extracts transcript sequences using gffread resulting in an
## "augmentad transcriptome" fasta file which is needed as a template for
## polyester in the next step.

set -e  # Exit on error

# ============================================================================
# CONFIGURATION
# ============================================================================

GTF_FILE="$1"
GENOME_FASTA="$2"
OUTPUT_FASTA="$3"

# ============================================================================
# MAIN EXECUTION
# ============================================================================

echo "=== Extracting Transcript Sequences ==="
echo ""

# Check if input files exist
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file not found: $GTF_FILE"
    exit 1
fi

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA file not found: $GENOME_FASTA"
    exit 1
fi

# Check if gffread is installed
if ! command -v gffread &> /dev/null; then
    echo "Error: gffread is not installed"
    echo "Install with: conda install -c bioconda gffread"
    exit 1
fi

# Extract transcript sequences
echo "Extracting transcript sequences with gffread..."
echo "  Input GTF: $GTF_FILE"
echo "  Genome: $GENOME_FASTA"
echo "  Output: $OUTPUT_FASTA"
echo ""

gffread -w "$OUTPUT_FASTA" -g "$GENOME_FASTA" "$GTF_FILE"

# Check output
if [ -f "$OUTPUT_FASTA" ]; then
    num_seqs=$(grep -c "^>" "$OUTPUT_FASTA")
    echo "Success! Extracted $num_seqs transcript sequences"
    echo "Output file: $OUTPUT_FASTA"
else
    echo "Error: Failed to create output file"
    exit 1
fi

echo ""
echo "=== Ready for Polyester simulation ==="
echo "Next step: Run polyester_cryptic_simulation.R"
