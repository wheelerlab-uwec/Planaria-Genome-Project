# assembly_stats_from_fasta.py
#
# Usage:
#   python assembly_stats_from_fasta.py assembly.fasta
#
# Reads a FASTA of contigs/scaffolds and reports:
#   - number of sequences
#   - total bp
#   - min / max / mean / median length
#   - N50, N90
#   - L50, L90

import sys
import statistics

def read_fasta_lengths(fasta_path):
    """
    Stream through a FASTA file and return a list of sequence lengths.
    Works even if sequences span multiple lines.
    """
    lengths = []
    current_len = 0

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # header line: commit previous sequence length (if any)
                if current_len > 0:
                    lengths.append(current_len)
                    current_len = 0
            else:
                # sequence line
                current_len += len(line)

    # last sequence
    if current_len > 0:
        lengths.append(current_len)

    return lengths

def compute_nx_and_lx(lengths, x_percent):
    """
    Compute Nx and Lx for a list of sequence lengths.

    Nx: contig length such that contigs of length >= Nx
        cover at least x_percent% of total assembly length.
    Lx: number of contigs contributing to that coverage.
    """
    if not lengths:
        return None, None

    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target = total_length * (x_percent / 100.0)

    cumulative = 0
    count = 0
    nx = None
    lx = None

    for length in sorted_lengths:
        cumulative += length
        count += 1
        if cumulative >= target:
            nx = length
            lx = count
            break

    return nx, lx

def format_bp(value):
    if value is None:
        return "NA"
    return f"{value:,} bp"

def main():
    if len(sys.argv) != 2:
        print("Usage: python assembly_stats_from_fasta.py assembly.fasta", file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]

    try:
        lengths = read_fasta_lengths(fasta_path)
    except FileNotFoundError:
        print(f"Error: FASTA file not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    if not lengths:
        print("Error: no sequences found in FASTA.", file=sys.stderr)
        sys.exit(1)

    num_seqs = len(lengths)
    total_bp = sum(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = total_bp / num_seqs
    median_len = statistics.median(lengths)

    n50, l50 = compute_nx_and_lx(lengths, 50)
    n90, l90 = compute_nx_and_lx(lengths, 90)

    print("Assembly statistics (from FASTA)")
    print("--------------------------------")
    print(f"FASTA file:      {fasta_path}")
    print()
    print(f"# sequences:     {num_seqs}")
    print(f"Total length:    {total_bp:,} bp")
    print(f"Min length:      {min_len:,} bp")
    print(f"Max length:      {max_len:,} bp")
    print(f"Mean length:     {mean_len:,.1f} bp")
    print(f"Median length:   {median_len:,.1f} bp")
    print()
    print(f"N50:             {format_bp(n50)}")
    print(f"L50:             {l50 if l50 is not None else 'NA'}")
    print(f"N90:             {format_bp(n90)}")
    print(f"L90:             {l90 if l90 is not None else 'NA'}")

if __name__ == "__main__":
    main()

