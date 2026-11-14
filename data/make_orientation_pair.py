"""Build the pUC19 orientation pair (3,686 bp each).

Downloads needed first (see README / Part 1 Section 3):
  - data/puc19.fasta: pUC19, GenBank L09137 (2,686 bp)
  - data/lambda.fasta: bacteriophage lambda, GenBank J02459 (48,502 bp)

  From NCBI Nucleotide, search the accession, then
  Send to -> File -> Format: FASTA -> Create File.

Method (in-silico model of a forward-vs-reverse cloning check):
  1. Find a lambda segment that has no BamHI site and exactly one
     site for another table enzyme, about one third from one end.
  2. Flank it with BamHI sites (GGATCC) so the inserted piece is
     1,000 bp in total (988 bp lambda + 2 x 6 bp BamHI ends).
  3. Insert it at the unique BamHI cut in pUC19, with BamHI sites
     regenerated on both sides as in a real ligation, once as
     written (forward) and once reverse-complemented (reverse).

Result: data/puc19_forward.fasta and data/puc19_reverse.fasta,
each 3,686 bp. BamHI alone gives 1,000 + 2,686 bp for both
(not distinguishing); pairing the insert's internal enzyme with
a backbone enzyme distinguishes them.

Usage (from the repository root):
    python data/make_orientation_pair.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from digest_planner.digest import find_sites
from digest_planner.enzymes import load_enzymes, reverse_complement
from digest_planner.sequence_io import Plasmid

INSERT_TOTAL = 1000
BAMHI_SITE = "GGATCC"
INNER_LEN = INSERT_TOTAL - 2 * len(BAMHI_SITE)  # 988 bp of lambda


def pick_lambda_segment(lambda_seq, enzymes):
    """Find a lambda window with no BamHI and one internal-enzyme site.

    Returns (start, internal_enzyme_name, site_start_in_insert).
    The internal site should sit about one third from one end so
    the forward/reverse fragments differ clearly.
    """
    candidates = []
    for start in range(0, len(lambda_seq) - INNER_LEN, 100):
        window = lambda_seq[start : start + INNER_LEN]
        if find_sites(window, BAMHI_SITE, circular=False):
            continue
        for name, enzyme in enzymes.items():
            if name == "BamHI":
                continue
            hits = find_sites(window, enzyme.site, circular=False)
            if len(hits) != 1:
                continue
            pos = hits[0]
            # About one third from one end (within 25%-42% of the window).
            frac = pos / INNER_LEN
            if 0.25 <= frac <= 0.42:
                candidates.append((start, name, pos))
                break
        if candidates:
            break
    if not candidates:
        raise RuntimeError(
            "No suitable lambda segment found: need no BamHI site and "
            "exactly one internal-enzyme site about a third from one end."
        )
    return candidates[0]


def main():
    root = Path(__file__).resolve().parent.parent
    data = root / "data"
    puc19_path = data / "puc19.fasta"
    lambda_path = data / "lambda.fasta"
    if not puc19_path.exists() or not lambda_path.exists():
        print("Missing input sequences. Download them first:")
        print("  - pUC19 (L09137) -> data/puc19.fasta")
        print("  - lambda (J02459) -> data/lambda.fasta")
        print("See README Section 'Datasets' for NCBI download steps.")
        return 2

    enzymes = load_enzymes(data / "enzymes.csv")
    puc19 = Plasmid.from_fasta(puc19_path)
    lam = Plasmid.from_fasta(lambda_path, circular=False)

    bamhi = enzymes["BamHI"]
    bam_sites = find_sites(puc19.sequence, bamhi.site, circular=True)
    if len(bam_sites) != 1:
        print(f"Expected exactly one BamHI site in pUC19, found {len(bam_sites)}.")
        return 2
    cut = (bam_sites[0] + bamhi.cut_offset) % puc19.length
    print(f"pUC19 length {puc19.length}; BamHI cut at {cut}")

    start, internal_name, internal_pos = pick_lambda_segment(lam.sequence, enzymes)
    inner = lam.sequence[start : start + INNER_LEN]
    insert_fwd = BAMHI_SITE + inner + BAMHI_SITE
    insert_rev = reverse_complement(insert_fwd)
    print(
        f"Lambda segment at {start} ({INNER_LEN} bp), "
        f"internal enzyme {internal_name} at offset {internal_pos}"
    )

    # Insert at the BamHI cut; BamHI sites are regenerated on both sides
    # because the insert itself carries BamHI ends.
    fwd = puc19.sequence[:cut] + insert_fwd + puc19.sequence[cut:]
    rev = puc19.sequence[:cut] + insert_rev + puc19.sequence[cut:]

    def write(path, name, seq):
        with open(path, "w") as handle:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i : i + 70] + "\n")

    write(data / "puc19_forward.fasta", "pUC19_forward insert fwd", fwd)
    write(data / "puc19_reverse.fasta", "pUC19_reverse insert rev", rev)
    print(f"Wrote {data / 'puc19_forward.fasta'} ({len(fwd)} bp)")
    print(f"Wrote {data / 'puc19_reverse.fasta'} ({len(rev)} bp)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
