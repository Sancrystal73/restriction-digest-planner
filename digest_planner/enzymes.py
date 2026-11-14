"""Restriction enzymes: parsing, validation, and loading.

Enzyme sites are written with '^' marking the top-strand cut, for
example EcoRI is G^AATTC. That means the site is GAATTC and the
cut comes after the first base (cut offset 1).

Only palindromic sites are kept: a site must equal its own
reverse complement.
"""

# Complement map: each base pairs with one other base.
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence.

    Example: reverse_complement("GGTCTC") gives "GAGACC".
    """
    result = ""
    for i in range(len(seq) - 1, -1, -1):
        result += COMPLEMENT[seq[i]]
    return result


def is_palindromic(site):
    """True when the site equals its own reverse complement."""
    return site == reverse_complement(site)


class Enzyme:
    """One restriction enzyme: a name, a site, and a cut offset."""

    def __init__(self, name, site, cut_offset):
        self.name = name
        self.site = site
        self.cut_offset = cut_offset

    def __str__(self):
        return self.name + " (" + self.site + ")"

    @classmethod
    def from_caret(cls, name, text):
        """Build an Enzyme from '^' notation, for example "CTGCA^G".

        The site is the text with '^' removed. The cut offset is
        the position of '^', the number of bases before the cut.
        """
        text = text.strip().upper()
        offset = text.find("^")
        if offset == -1:
            raise ValueError("Enzyme '" + name + "': site '" + text
                             + "' has no '^' marking the cut")
        site = text.replace("^", "")
        return cls(name, site, offset)


def load_enzymes(path):
    """Load an enzyme table from a CSV file.

    The file has a header line (name,site) and one row per enzyme,
    with the site in '^' notation. Rows with a non-palindromic
    site are skipped. Returns a dictionary of name -> Enzyme.
    """
    enzymes = {}
    with open(path) as f:
        f.readline()  # skip the header line
        for line in f:
            line = line.strip()
            if line == "":
                continue
            name, site_text = line.split(",")
            enzyme = Enzyme.from_caret(name.strip(), site_text.strip())
            if not is_palindromic(enzyme.site):
                print("Warning: skipping non-palindromic enzyme "
                      + name + " (site " + enzyme.site + ")")
                continue
            enzymes[name] = enzyme
    return enzymes
