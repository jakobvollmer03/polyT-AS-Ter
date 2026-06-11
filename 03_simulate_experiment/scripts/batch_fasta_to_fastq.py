import os
import gzip


def convert(fasta, fastq):
    """Convert FASTA to FASTQ with uniform quality scores."""
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {fasta}")
    
    with open(fasta) as f_in, gzip.open(fastq, "wt") as f_out:
        header = None
        seq_lines = []

        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seq = "".join(seq_lines)
                    qual = "I" * len(seq)
                    f_out.write(f"@{header}\n{seq}\n+\n{qual}\n")
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            seq = "".join(seq_lines)
            qual = "I" * len(seq)
            f_out.write(f"@{header}\n{seq}\n+\n{qual}\n")


# Get values from Snakemake
fasta_dir = snakemake.input[0]  # Directory path from input
sample = snakemake.wildcards.sample
prefix = snakemake.params.prefix

# Construct input fasta paths
fasta1 = os.path.join(fasta_dir, f"{prefix}_sample_{sample}_1.fasta")
fasta2 = os.path.join(fasta_dir, f"{prefix}_sample_{sample}_2.fasta")

# Output paths already defined by Snakemake
fastq1 = snakemake.output[0]  # fastq1
fastq2 = snakemake.output[1]  # fastq2

# Ensure output directory exists
os.makedirs(os.path.dirname(fastq1), exist_ok=True)

# Convert
convert(fasta1, fastq1)
convert(fasta2, fastq2)