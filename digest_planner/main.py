"""Command-line entry point for the Restriction Digest Planner.

Usage:
    python -m digest_planner.main plasmid_A.fasta plasmid_B.fasta enzymes.csv
    python -m digest_planner.main plasmid_A.fasta plasmid_B.fasta enzymes.csv 0.10 5

The first three arguments are the two FASTA files and the enzyme
table. The optional fourth argument is the version 2 threshold
(for example 0.10); leaving it out runs version 1, which reports
every digest whose fragment lists differ at all. The optional
fifth argument is how many top digests to report (default 5).
"""

import sys

from .enzymes import load_enzymes
from .planner import DigestPlanner
from .sequence_io import Plasmid


def main():
    args = sys.argv[1:]
    if len(args) < 3 or len(args) > 5:
        print("Usage: python -m digest_planner.main <plasmid_A.fasta> "
              "<plasmid_B.fasta> <enzymes.csv> [threshold] [top]")
        return 2
    threshold = None
    if len(args) >= 4:
        threshold = float(args[3])
        if threshold < 0.0 or threshold > 1.0:
            print("Threshold must be between 0 and 1")
            return 2
    top = 5
    if len(args) >= 5:
        top = int(args[4])
    plasmid_a = Plasmid.from_fasta(args[0])
    plasmid_b = Plasmid.from_fasta(args[1])
    enzymes = load_enzymes(args[2])
    if len(enzymes) == 0:
        print("No usable enzymes loaded from " + args[2])
        return 2
    planner = DigestPlanner(plasmid_a, plasmid_b, enzymes, threshold)
    planner.run()
    planner.report(top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
