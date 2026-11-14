"""Generate the synthetic test pair (3,000 bp each).

Filler bases are drawn at random from a fixed seed. Any filler
sequence containing a site from the enzyme table is cleaned by
mutating one base inside each unintended site (away from the
planted sites), so only the planted sites exist in the end.

Planted sites reproduce the worked example in Part 1:
  - EcoRI (GAATTC, cut offset 1) cuts at position 100 in both
    plasmids, so its site starts at 99.
  - BamHI (GGATCC, cut offset 1) cuts at position 900 in A
    (site start 899) and at position 1,500 in B (site start 1499).

Expected digests:
  - EcoRI alone / BamHI alone: [3000] vs [3000], not distinguishing.
  - EcoRI + BamHI: [800, 2200] vs [1400, 1600], distinguishing,
    score 0.43.

Usage (from the repository root):
    python data/make_synthetic.py
    python data/make_synthetic.py --seed 42 --outdir data
"""

import argparse
import random
import sys
from pathlib import Path

# Allow running as `python data/make_synthetic.py` from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from digest_planner.digest import find_sites
from digest_planner.enzymes import load_enzymes

LENGTH = 3000
BASES = "ACGT"

# (enzyme name, site, site start) for each plasmid.
PLANTED_A = [("EcoRI", "GAATTC", 99), ("BamHI", "GGATCC", 899)]
PLANTED_B = [("EcoRI", "GAATTC", 99), ("BamHI", "GGATCC", 1499)]


def _plant(seq_list, site, start):
    for i, base in enumerate(site):
        seq_list[(start + i) % len(seq_list)] = base


def _planted_ranges(planted, length):
    """Set of positions covered by planted sites (for protection)."""
    covered = set()
    for _, site, start in planted:
        for i in range(len(site)):
            covered.add((start + i) % length)
    return covered


def _clean_unintended(seq_list, enzymes, planted, rng, max_rounds=500):
    """Mutate unintended sites until only planted sites remain."""
    length = len(seq_list)
    protected = _planted_ranges(planted, length)
    planted_set = {(site, start) for _, site, start in planted}
    for _ in range(max_rounds):
        seq = "".join(seq_list)
        bad = []
        for enzyme in enzymes.values():
            for start in find_sites(seq, enzyme.site, circular=True):
                if (enzyme.site, start) not in planted_set:
                    bad.append((enzyme.site, start))
        if not bad:
            return seq_list
        # Mutate one base in the first unintended site found, picking
        # a position that is not part of a planted site.
        site, start = bad[0]
        candidates = [
            (start + i) % length
            for i in range(len(site))
            if (start + i) % length not in protected
        ]
        if not candidates:
            # Unintended site overlaps a planted site; nudge a
            # neighbouring unprotected base to break it.
            candidates = [(start + len(site)) % length]
        pos = candidates[len(candidates) // 2]
        current = seq_list[pos]
        choices = [b for b in BASES if b != current]
        seq_list[pos] = rng.choice(choices)
    raise RuntimeError("Could not clean unintended sites after many rounds")


def make_plasmid(planted, enzymes, rng):
    seq_list = [rng.choice(BASES) for _ in range(LENGTH)]
    for _, site, start in planted:
        _plant(seq_list, site, start)
    _clean_unintended(seq_list, enzymes, planted, rng)
    # Re-plant (cleaning never touches protected positions, but be safe).
    for _, site, start in planted:
        _plant(seq_list, site, start)
    return "".join(seq_list)


def write_fasta(path, name, seq, width=70):
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(seq), width):
            handle.write(seq[i : i + width] + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate synthetic plasmid pair.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default=str(Path(__file__).parent))
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parent.parent
    enzymes = load_enzymes(root / "data" / "enzymes.csv")
    rng = random.Random(args.seed)

    seq_a = make_plasmid(PLANTED_A, enzymes, rng)
    seq_b = make_plasmid(PLANTED_B, enzymes, rng)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_fasta(outdir / "plasmid_A.fasta", "plasmid_A synthetic 3000bp", seq_a)
    write_fasta(outdir / "plasmid_B.fasta", "plasmid_B synthetic 3000bp", seq_b)
    print(f"Wrote {outdir / 'plasmid_A.fasta'} ({len(seq_a)} bp)")
    print(f"Wrote {outdir / 'plasmid_B.fasta'} ({len(seq_b)} bp)")

    # Verify: only the planted sites exist.
    for label, seq, planted in [("A", seq_a, PLANTED_A), ("B", seq_b, PLANTED_B)]:
        found = []
        for enzyme in enzymes.values():
            for start in find_sites(seq, enzyme.site, circular=True):
                found.append((enzyme.name, start))
        print(f"Plasmid {label} sites: {sorted(found)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
