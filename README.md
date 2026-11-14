# Restriction Digest Planner — CS 313E Term Project, Part 2 (Implementation)

**Choosing enzymes to distinguish two candidate plasmids**

After a cloning step, a researcher often holds a plasmid that could be one
of two constructs — for example the same insert ligated in the forward or
in the reverse orientation. Both candidates have the same length, so an
uncut or once-cut plasmid looks identical on an agarose gel. This program
finds which single- or double-enzyme digest produces gel patterns that
tell the two candidates apart.

Given two plasmid sequences (FASTA) and a table of available enzymes
(name, recognition sequence, cut position), the program:

1. Locates every recognition site in each plasmid, treating each plasmid
   as a circle (sites may span the numbering origin).
2. Simulates every single-enzyme digest and every two-enzyme digest and
   computes the expected fragment lengths for plasmid A and plasmid B.
3. Reports which digests distinguish A from B, the expected fragments
   for each, and a ranking by how easily the bands could be told apart.

Part 1 (project definition) is in the course submission; this repository
is the Part 2 implementation. It is written in Python 3 using only the
standard library.

## Quick start

```bash
# Generate the synthetic test pair (3,000 bp each)
python data/make_synthetic.py

# Version 1: report every digest whose fragment lists differ at all
python -m digest_planner.main data/plasmid_A.fasta data/plasmid_B.fasta data/enzymes.csv --top 5

# Version 2 (improved): only digests with a >=10% band difference
python -m digest_planner.main data/plasmid_A.fasta data/plasmid_B.fasta data/enzymes.csv --threshold 0.10 --top 5

# Run the tests (34 tests; 2 skipped unless pUC19/lambda are downloaded)
python -m unittest discover -s tests
```

No third-party packages are needed to run the program. Python 3.10+
is sufficient.

## Example run

Synthetic pair, 12 enzymes, 78 digests tried:

```
$ python -m digest_planner.main data/plasmid_A.fasta data/plasmid_B.fasta data/enzymes.csv --threshold 0.10 --top 5

Plasmid A: plasmid_A (3000 bp)
Plasmid B: plasmid_B (3000 bp)
Enzymes: 12; digests tried: 78; distinguishing (v2 (threshold=0.1)): 1

Digest                    Plasmid A (bp)       Plasmid B (bp)       Score
BamHI + EcoRI             800, 2200            1400, 1600           0.43
```

Saved outputs are in `results/`.

## Results and interpretation

### Synthetic pair (3,000 bp each)

| Digest | Plasmid A (bp) | Plasmid B (bp) | Distinguishable? |
| --- | --- | --- | --- |
| EcoRI only | 3,000 | 3,000 | No |
| BamHI only | 3,000 | 3,000 | No |
| EcoRI + BamHI | 800, 2,200 | 1,400, 1,600 | **Yes (score 0.43)** |

Each enzyme cuts each plasmid once, so either single digest only opens
the circle into one 3,000 bp piece — the two candidates look identical.
The double digest measures the distance between the two sites, which is
800 bp in A but 1,400 bp in B, so only the double digest reveals the
difference. That is exactly the worked example from Part 1, and the
program reproduces it: of all 78 digests, only BamHI + EcoRI
distinguishes the pair. Every other enzyme cuts neither plasmid (no
sites by construction), and every other pair involves at most one
cutting enzyme, so it behaves like a single digest.

The full run finishes in about 0.08 s on a laptop, well under the
one-second expectation from Part 1.

### What the score means

The score is the largest nearest-band relative difference between the
two lanes. A score of 0.43 for BamHI + EcoRI means at least one band
differs by 43% from every band in the other lane (the 800 bp band in A
vs the 1,400 bp band in B) — an easy call on a real gel. A score near
zero means the lanes are effectively identical.

### Improvement applied: version 2 threshold (v1 → v2)

Version 1 reports a digest as distinguishing whenever the fragment
lists differ at all. That over-reports: on a real gel, bands a few
percent apart run as one band, because migration distance is roughly
linear in the *logarithm* of fragment size.

Version 2 adds an adjustable relative-difference threshold
(`--threshold`, default 0.10). A digest counts as distinguishing only
when its score reaches the threshold, i.e. at least one band differs
by that fraction from every band in the other lane.

Concrete case (test T12): fragments `[1000, 2000]` vs `[1020, 1980]`.

- Version 1: distinguishing (the lists differ).
- Version 2 at threshold 0.05: **not** distinguishing (score 0.0196 —
  the bands differ by only ~2%, which would appear as single bands).

On the synthetic pair, v1 and v2 agree (the one distinguishing digest
scores 0.43, far above any of 0.05 / 0.10 / 0.15), which is the expected
outcome for a clean, well-separated example. The improvement matters on
realistic inputs, where near-identical fragment lists are common: v2
reports fewer digests, but ones more likely to be readable on a gel.
Part 3 will compare the two versions on both datasets at thresholds
0.05, 0.10, and 0.15.

Other improvements applied in this implementation:

- **One-uncut handling:** if exactly one plasmid is uncut, the result
  is distinguishing but scored 0.00 and ranked last, because uncut
  (supercoiled / nicked) plasmid does not run at its true size.
- **Duplicate-cut removal:** two enzymes sharing a site and cut
  position produce one cut, never a 0 bp fragment.
- **Cut-position caching:** each single-enzyme search runs once per
  plasmid and is reused by all 11 pair digests involving that enzyme.

### Expected weaknesses (from Part 1, still true)

- Only palindromic sites written with A, C, G, T are supported;
  ambiguous bases and Type IIS enzymes are excluded.
- At most two enzymes per digest.
- Biochemistry is not modelled: Dam/Dcm methylation, incomplete
  digestion, and double-digest buffer compatibility are all assumed away.
- A single relative threshold stands in for gel resolution; band
  intensity, very small fragments, and the 100 bp–25 kb usable range
  are not modelled. The 0.10 default is a working assumption, not a
  measured value.
- Selection sort (O(n²)) and brute-force search (O(L·m)) suit plasmids
  but would be slow on genome-sized input.

## Datasets

| Dataset | Contents | Source |
| --- | --- | --- |
| Synthetic test pair | Plasmids A and B, 3,000 bp each | Generated by `data/make_synthetic.py` (fixed seed 42); included as `data/plasmid_A.fasta`, `data/plasmid_B.fasta` |
| pUC19 orientation pair | pUC19 + 1,000 bp lambda insert, forward (A) and reverse (B), 3,686 bp each | Built by `data/make_orientation_pair.py` from the two public sequences below |
| pUC19 cloning vector | 2,686 bp | GenBank L09137 |
| Bacteriophage lambda | 48,502 bp (source of the 1,000 bp insert) | GenBank J02459 |
| Enzyme table | 12 common Type IIP enzymes | `data/enzymes.csv`, compiled by hand from REBASE |

**Downloading the public sequences** (needed only for the orientation
pair and tests T14/T16):

1. Open [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) and
   search for the accession (L09137 or J02459).
2. Choose Send to → File → Format: FASTA → Create File.
3. Save the files as `data/puc19.fasta` and `data/lambda.fasta`.
4. Run `python data/make_orientation_pair.py`.

**Enzyme table format.** One row per enzyme; `^` marks the top-strand cut:

```csv
name,site
EcoRI,G^AATTC
BamHI,G^GATCC
```

## Design

Six modules in a pipeline: input parsing → digest simulation →
comparison → ranking.

```
digest_planner/
  main.py         command-line entry point
  sequence_io.py  FASTA reading and validation
  enzymes.py      Enzyme class, CSV loading, reverse complement
  digest.py       site search, selection sort, fragment lengths
  compare.py      distinguishability tests and scoring
  planner.py      DigestPlanner: enumerate, evaluate, rank, report
data/
  enzymes.csv, plasmid_A.fasta, plasmid_B.fasta,
  make_synthetic.py, make_orientation_pair.py
tests/
  test_enzymes.py, test_digest.py, test_compare.py, test_planner.py
results/
  example outputs for the synthetic pair
```

| Module | Key classes / functions | Responsibility |
| --- | --- | --- |
| `sequence_io.py` | `Plasmid`, `Plasmid.from_fasta`, `validate` | Reads one FASTA record, upper-cases it, rejects non-ACGT |
| `enzymes.py` | `Enzyme`, `Enzyme.from_caret`, `is_palindromic`, `reverse_complement`, `load_enzymes` | Parses the CSV, converts `^` notation to site + offset, keeps palindromic enzymes |
| `digest.py` | `find_sites`, `selection_sort`, `cut_positions`, `fragment_lengths` | Wrap-around site search, sorted/deduped cuts, fragments |
| `compare.py` | `differs_exact`, `separation_score`, `distinguishable` | Decides whether two fragment lists look different; scores them |
| `planner.py` | `DigestResult`, `DigestPlanner` (`enumerate_schemes`, `evaluate`, `run`, `top`, `report`) | Searches all digests, ranks the distinguishing ones, prints the report |
| `main.py` | `main` | Parses arguments, runs the planner |

Pipeline: `main` loads the plasmids and enzymes → `enumerate_schemes`
lists every single enzyme and unordered pair (78 for 12 enzymes) →
for each scheme and plasmid, `cut_positions` finds sites, converts
them to cut positions, sorts and dedupes them, and `fragment_lengths`
turns them into fragments → `evaluate` compares the two fragment
lists → `run` pushes distinguishing results onto a heap keyed by
score → `top(k)` pops the k best and `report` prints them.

### Fragments from a circular plasmid

With cut positions sorted as c₁ < c₂ < … < cₙ in a plasmid of length L:

- 0 cuts: uncut; 1 cut: one fragment of length L.
- Otherwise fᵢ = cᵢ₊₁ − cᵢ, and fₙ = (L − cₙ) + c₁ (wraps the origin).
  All lengths sum to L.

### Libraries

Standard library only: `argparse` (CLI), `csv` (enzyme table), `heapq`
(ranking), `random` (synthetic data, fixed seed), `unittest` (tests).
Biopython (`Bio.Restriction`) is used only in the optional cross-check
test T16 and is not required to run the program.

### Data structures and algorithms

- **Dictionary:** enzyme table (name → Enzyme), complement map,
  cut-position cache per enzyme and plasmid.
- **List:** site positions, cut positions, fragment lengths, schemes.
- **Heap (priority queue):** ranking distinguishing digests by score.
- **Selection sort:** ordering cut positions and fragment lists.
- **Brute-force string search:** sliding-window site search, extended
  by m−1 bases for circular wrap-around.
- **Exhaustive enumeration:** all singles and unordered pairs.
- **Heap-based selection:** popping the k highest scores for the report.

## Tests

16 tests from Part 1 (T1–T16), written with `unittest`:

- T1–T3: reverse complement, palindrome filter, `^` parsing.
- T4–T10: site search (linear and origin-spanning), selection sort,
  fragment lengths, zero/one cut, duplicate cuts, fragment-sum invariant
  (100 random sequences).
- T11–T12: version 1 rule and version 2 threshold.
- T13: end-to-end synthetic pair (single digests fail, double succeeds,
  score 0.43).
- T14: end-to-end orientation pair (skipped unless the pair is generated).
- T15: FASTA validation rejects non-ACGT characters.
- T16: optional Biopython cross-check (skipped unless pUC19 is
  downloaded and Biopython is installed).

```bash
python -m unittest discover -s tests
# Ran 34 tests in ~0.5s — OK (skipped=2)
```

## References

1. Addgene. *Diagnostic Restriction Digest* (protocol).
   <https://www.addgene.org/protocols/diagnostic-digest/>
2. Addgene. *Restriction Digest of Plasmid DNA* (protocol).
   <https://www.addgene.org/protocols/restriction-digest/>
3. Lee, P. Y., Costumbrado, J., Hsu, C.-Y., and Kim, Y. H. (2012).
   Agarose gel electrophoresis for the separation of DNA fragments.
   *Journal of Visualized Experiments*, (62), 3923. doi:10.3791/3923.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC4846332>
4. O'Connor, C. *Investigations in Molecular Cell Biology*, Section 11.1:
   Restriction endonucleases. LibreTexts Biology.
5. Roberts, R. J., Vincze, T., Posfai, J., and Macelis, D. (2010). REBASE.
   *Nucleic Acids Research*, 38 (Database issue), D234–D236.
   Database: <http://rebase.neb.com>
6. SnapGene. *Restriction Enzymes: A Brief History* (guide).
   <https://www.snapgene.com/guides/restriction-enzymes>
7. New England Biolabs. *pUC19 map and sequence* (GenBank L09137).
8. European Nucleotide Archive. *Escherichia phage Lambda, complete
   genome* (J02459). <https://www.ebi.ac.uk/ena/browser/view/J02459>
9. National Center for Biotechnology Information. *Nucleotide database*.
   <https://www.ncbi.nlm.nih.gov/nucleotide>
