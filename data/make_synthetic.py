"""Generate the synthetic test pair (3,000 bp each).

Filler bases are drawn at random from a fixed seed. Any filler
sequence containing a site from the enzyme table is cleaned by
changing one base inside each unintended site (away from the
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
"""

import os
import random
import sys

# Let this script import the package when run from the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from digest_planner.digest import find_sites
from digest_planner.enzymes import load_enzymes

LENGTH = 3000
BASES = "ACGT"

# (enzyme name, site, site start) for each plasmid.
PLANTED_A = [("EcoRI", "GAATTC", 99), ("BamHI", "GGATCC", 899)]
PLANTED_B = [("EcoRI", "GAATTC", 99), ("BamHI", "GGATCC", 1499)]


def plant(seq, site, start):
    """Write a site into the sequence list at a start position."""
    for i in range(len(site)):
        seq[(start + i) % len(seq)] = site[i]


def clean_unintended(seq, enzymes, planted, rng):
    """Change bases until only the planted sites remain."""
    length = len(seq)
    protected = set()
    for name, site, start in planted:
        for i in range(len(site)):
            protected.add((start + i) % length)
    planted_set = set()
    for name, site, start in planted:
        planted_set.add((site, start))
    for _ in range(500):
        text = "".join(seq)
        bad = None
        for enzyme in enzymes.values():
            for start in find_sites(text, enzyme.site, True):
                if (enzyme.site, start) not in planted_set:
                    bad = (enzyme.site, start)
                    break
            if bad is not None:
                break
        if bad is None:
            return seq
        site, start = bad
        # Change a base inside the unintended site that is not
        # part of a planted site.
        pos = None
        for i in range(len(site)):
            p = (start + i) % length
            if p not in protected:
                pos = p
                break
        if pos is None:
            pos = (start + len(site)) % length
        choices = [b for b in BASES if b != seq[pos]]
        seq[pos] = rng.choice(choices)
    raise RuntimeError("Could not clean unintended sites")


def make_plasmid(planted, enzymes, rng):
    """Build one 3,000 bp plasmid with only the planted sites."""
    seq = [rng.choice(BASES) for _ in range(LENGTH)]
    for name, site, start in planted:
        plant(seq, site, start)
    clean_unintended(seq, enzymes, planted, rng)
    for name, site, start in planted:
        plant(seq, site, start)
    return "".join(seq)


def write_fasta(path, name, seq):
    """Write a sequence to a FASTA file, 70 bases per line."""
    with open(path, "w") as f:
        f.write(">" + name + "\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    enzymes = load_enzymes(os.path.join(root, "data", "enzymes.csv"))
    rng = random.Random(42)
    seq_a = make_plasmid(PLANTED_A, enzymes, rng)
    seq_b = make_plasmid(PLANTED_B, enzymes, rng)
    data_dir = os.path.join(root, "data")
    write_fasta(os.path.join(data_dir, "plasmid_A.fasta"),
                "plasmid_A synthetic 3000bp", seq_a)
    write_fasta(os.path.join(data_dir, "plasmid_B.fasta"),
                "plasmid_B synthetic 3000bp", seq_b)
    print("Wrote plasmid_A.fasta (" + str(len(seq_a)) + " bp)")
    print("Wrote plasmid_B.fasta (" + str(len(seq_b)) + " bp)")
    for label, seq in [("A", seq_a), ("B", seq_b)]:
        found = []
        for enzyme in enzymes.values():
            for start in find_sites(seq, enzyme.site, True):
                found.append((enzyme.name, start))
        print("Plasmid " + label + " sites: " + str(sorted(found)))


if __name__ == "__main__":
    main()
