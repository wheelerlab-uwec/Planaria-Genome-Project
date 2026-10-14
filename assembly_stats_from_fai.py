# assembly_stats_from_fai.py
#
# Usage:
#   python assembly_stats_from_fai.py draft_assembly.fasta.fai
#
# Reads a samtools-style .fai index and reports:
#   - number of sequences
#   - total bp
#   - min / max / mean / median length
#   - N50, N90
#   - L50, L90
#
# This avoids loading the full FASTA into memory by using just the index.

import sys
import statistics

def compute_nx_and_lx(lengths, x_percent):
    """
    Compute Nx and Lx for a list of sequence lengths.

    Nx: contig length such that contigs of length >= Nx
        cover at least x_percent% of the total assembly length.
    Lx: number of contigs contributing to that coverage.

    Returns (Nx, Lx). If lengths is empty, returns (None, None).
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


def load_lengths_from_fai(fai_path):
    """
    Parse a .fai file and return a list of sequence lengths.
    Assumes standard samtools FAI format:
       <name>\t<length>\t<offset>\t<bases_per_line>\t<bytes_per_line>
    """
    lengths = []

    with open(fai_path, "r") as fai_file:
        for line in fai_file:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                length = int(parts[1])
                lengths.append(length)
            except ValueError:
                # Skip malformed lines rather than failing hard
                continue

    return lengths


def format_bp(value):
    """Return an integer bp value with thousands separators, or 'NA'."""
    if value is None:
        return "NA"
    return f"{value:,} bp"


def main():
    if len(sys.argv) != 2:
        print("Usage: python assembly_stats_from_fai.py draft_assembly.fasta.fai", file=sys.stderr)
        sys.exit(1)

    fai_path = sys.argv[1]

    try:
        lengths = load_lengths_from_fai(fai_path)
    except FileNotFoundError:
        print(f"Error: FAI file not found: {fai_path}", file=sys.stderr)
        sys.exit(1)

    if not lengths:
        print("Error: no valid sequence lengths found in FAI file.", file=sys.stderr)
        sys.exit(1)

    # Basic stats
    num_seqs = len(lengths)
    total_bp = sum(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = total_bp / num_seqs
    median_len = statistics.median(lengths)

    # Nx/Lx stats
    n50, l50 = compute_nx_and_lx(lengths, 50)
    n90, l90 = compute_nx_and_lx(lengths, 90)

    print("Assembly statistics (from FAI)")
    print("--------------------------------")
    print(f"FAI file:        {fai_path}")
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

