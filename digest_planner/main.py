"""Command-line entry point for the Restriction Digest Planner.

Usage:
    python -m digest_planner.main data/plasmid_A.fasta data/plasmid_B.fasta \\
        data/enzymes.csv --threshold 0.10 --top 5

Omitting --threshold runs version 1, which reports every digest
whose fragment lists differ at all. With --threshold, version 2
reports only digests whose separation score reaches the threshold.
"""

import argparse
import sys

# Support both `python -m digest_planner.main` and `python main.py`.
try:
    from .enzymes import load_enzymes
    from .planner import DigestPlanner
    from .sequence_io import Plasmid
except ImportError:  # pragma: no cover - direct script run
    from enzymes import load_enzymes
    from planner import DigestPlanner
    from sequence_io import Plasmid


def build_parser():
    parser = argparse.ArgumentParser(
        description="Choose restriction digests that distinguish two plasmids."
    )
    parser.add_argument("plasmid_a", help="FASTA file for plasmid A")
    parser.add_argument("plasmid_b", help="FASTA file for plasmid B")
    parser.add_argument("enzymes", help="Enzyme table CSV (name,site with ^ cut mark)")
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Version 2 relative-difference threshold, e.g. 0.10. "
        "Omit for version 1 (any difference counts).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="How many top distinguishing digests to report (default 5).",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.threshold is not None and not (0.0 <= args.threshold <= 1.0):
        print("--threshold must be between 0 and 1", file=sys.stderr)
        return 2
    plasmid_a = Plasmid.from_fasta(args.plasmid_a)
    plasmid_b = Plasmid.from_fasta(args.plasmid_b)
    enzymes = load_enzymes(args.enzymes)
    if not enzymes:
        print(f"No usable enzymes loaded from {args.enzymes}", file=sys.stderr)
        return 2
    planner = DigestPlanner(plasmid_a, plasmid_b, enzymes, threshold=args.threshold)
    planner.run()
    planner.report(k=args.top)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
