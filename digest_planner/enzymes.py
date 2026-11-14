"""Restriction enzymes: parsing, validation, and loading.

Enzyme sites are written with '^' marking the top-strand cut, e.g.
EcoRI is ``G^AATTC``. The loader converts that into a recognition
site (``GAATTC``) and a cut offset (1, the number of bases before
the cut). Only palindromic sites are kept: a site must equal its
own reverse complement.
"""

import csv
import warnings
from pathlib import Path


# Complement map (dictionary data structure from the course).
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence.

    Example: reverse_complement("GGTCTC") -> "GAGACC".
    """
    seq = seq.upper()
    return "".join(_COMPLEMENT[base] for base in reversed(seq))


def is_palindromic(site):
    """True when the site equals its own reverse complement."""
    site = site.upper()
    return site == reverse_complement(site)


class Enzyme:
    """A Type IIP restriction enzyme.

    Attributes:
        name: Enzyme name, e.g. "EcoRI".
        site: Recognition sequence, e.g. "GAATTC".
        cut_offset: Number of bases before the top-strand cut.
            For EcoRI (G^AATTC) this is 1.
    """

    def __init__(self, name, site, cut_offset):
        self.name = name
        self.site = site.upper()
        self.cut_offset = int(cut_offset)

    @classmethod
    def from_caret(cls, name, text):
        """Build an Enzyme from '^' notation, e.g. "CTGCA^G".

        The site is the text with '^' removed; the cut offset is
        the index of '^'.
        """
        text = text.strip().upper()
        if "^" not in text:
            raise ValueError(
                f"Enzyme {name!r}: site {text!r} has no '^' marking the cut"
            )
        offset = text.index("^")
        site = text.replace("^", "")
        return cls(name, site, offset)

    def __repr__(self):
        return f"Enzyme(name={self.name!r}, site={self.site!r}, cut_offset={self.cut_offset})"


def load_enzymes(path):
    """Load an enzyme table from a CSV file.

    The CSV has a header row (name, site) and one row per enzyme,
    with the site in '^' notation. Non-palindromic rows are
    skipped with a warning. Returns a dict of name -> Enzyme.
    """
    enzymes = {}
    path = Path(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = (row.get("name") or "").strip()
            site_text = (row.get("site") or "").strip()
            if not name or not site_text:
                continue
            enzyme = Enzyme.from_caret(name, site_text)
            if not is_palindromic(enzyme.site):
                warnings.warn(
                    f"Skipping non-palindromic enzyme {name}: site {enzyme.site}",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            enzymes[name] = enzyme
    return enzymes
