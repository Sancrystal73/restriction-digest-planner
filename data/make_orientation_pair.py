"""Build the pUC19 orientation pair (3,686 bp each).

Downloads needed first (see README):
  - data/puc19.fasta: pUC19, GenBank L09137 (2,686 bp)
  - data/lambda.fasta: bacteriophage lambda, GenBank J02459 (48,502 bp)

  From NCBI Nucleotide, search the accession, then
  Send to -> File -> Format: FASTA -> Create File.

Method (a model of a forward-vs-reverse cloning check):
  1. Find a lambda segment with no BamHI site and exactly one
     site for another table enzyme, about one third from one end.
  2. Put BamHI sites (GGATCC) on both ends so the inserted piece
     is 1,000 bp in total (988 bp lambda + 2 x 6 bp ends).
  3. Insert it at the unique BamHI cut in pUC19, once as written
     (forward) and once reverse-complemented (reverse).

Result: data/puc19_forward.fasta and data/puc19_reverse.fasta,
each 3,686 bp. BamHI alone gives 1,000 + 2,686 bp for both
(not distinguishing); pairing the insert's internal enzyme with
a backbone enzyme distinguishes them.

Usage (from the repository root):
    python data/make_orientation_pair.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from digest_planner.digest import find_sites
from digest_planner.enzymes import load_enzymes, reverse_complement
from digest_planner.sequence_io import Plasmid

INSERT_TOTAL = 1000
BAMHI_SITE = "GGATCC"
INNER_LEN = INSERT_TOTAL - 2 * len(BAMHI_SITE)  # 988 bp of lambda


def pick_lambda_segment(lambda_seq, enzymes):
    """Find a lambda window with no BamHI and one internal site.

    Returns (start, enzyme name, site position in the window).
    The internal site should sit about one third from one end so
    the forward/reverse fragments differ clearly.
    """
    for start in range(0, len(lambda_seq) - INNER_LEN, 100):
        window = lambda_seq[start:start + INNER_LEN]
        if len(find_sites(window, BAMHI_SITE, False)) > 0:
            continue
        for name in enzymes:
            if name == "BamHI":
                continue
            enzyme = enzymes[name]
            hits = find_sites(window, enzyme.site, False)
            if len(hits) != 1:
                continue
            frac = hits[0] / INNER_LEN
            if 0.25 <= frac <= 0.42:
                return (start, name, hits[0])
    raise RuntimeError("No suitable lambda segment found")


def write_fasta(path, name, seq):
    """Write a sequence to a FASTA file, 70 bases per line."""
    with open(path, "w") as f:
        f.write(">" + name + "\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(root, "data")
    puc19_path = os.path.join(data_dir, "puc19.fasta")
    lambda_path = os.path.join(data_dir, "lambda.fasta")
    if not os.path.exists(puc19_path) or not os.path.exists(lambda_path):
        print("Missing input sequences. Download them first:")
        print("  - pUC19 (L09137) -> data/puc19.fasta")
        print("  - lambda (J02459) -> data/lambda.fasta")
        print("See README for the NCBI download steps.")
        return 2

    enzymes = load_enzymes(os.path.join(data_dir, "enzymes.csv"))
    puc19 = Plasmid.from_fasta(puc19_path)
    lam = Plasmid.from_fasta(lambda_path, False)

    bamhi = enzymes["BamHI"]
    bam_sites = find_sites(puc19.sequence, bamhi.site, True)
    if len(bam_sites) != 1:
        print("Expected exactly one BamHI site in pUC19, found "
              + str(len(bam_sites)) + ".")
        return 2
    cut = (bam_sites[0] + bamhi.cut_offset) % puc19.length
    print("pUC19 length " + str(puc19.length)
          + "; BamHI cut at " + str(cut))

    start, internal_name, internal_pos = pick_lambda_segment(lam.sequence,
                                                             enzymes)
    inner = lam.sequence[start:start + INNER_LEN]
    insert_fwd = BAMHI_SITE + inner + BAMHI_SITE
    insert_rev = reverse_complement(insert_fwd)
    print("Lambda segment at " + str(start) + " (" + str(INNER_LEN)
          + " bp), internal enzyme " + internal_name
          + " at offset " + str(internal_pos))

    # Insert at the BamHI cut. BamHI sites are regenerated on both
    # sides because the insert itself carries BamHI ends.
    fwd = puc19.sequence[:cut] + insert_fwd + puc19.sequence[cut:]
    rev = puc19.sequence[:cut] + insert_rev + puc19.sequence[cut:]
    write_fasta(os.path.join(data_dir, "puc19_forward.fasta"),
                "pUC19_forward insert fwd", fwd)
    write_fasta(os.path.join(data_dir, "puc19_reverse.fasta"),
                "pUC19_reverse insert rev", rev)
    print("Wrote puc19_forward.fasta (" + str(len(fwd)) + " bp)")
    print("Wrote puc19_reverse.fasta (" + str(len(rev)) + " bp)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
