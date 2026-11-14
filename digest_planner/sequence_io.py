"""FASTA reading and validation for plasmid sequences.

A plasmid is stored as a Plasmid object with a name, an uppercase
ACGT-only sequence, its length, and whether it is circular.
"""

from pathlib import Path


VALID_BASES = set("ACGT")


class Plasmid:
    """A single plasmid sequence.

    Attributes:
        name: Identifier taken from the FASTA header (first word).
        sequence: Uppercase DNA sequence, ACGT only.
        length: Length of the sequence in bp.
        circular: True for plasmids (default); False for linear DNA.
    """

    def __init__(self, name, sequence, circular=True):
        self.name = name
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.circular = circular
        self.validate()

    def validate(self):
        """Reject any character other than A, C, G, T.

        Raises:
            ValueError: naming the offending character and its
                0-based position in the sequence.
        """
        for i, base in enumerate(self.sequence):
            if base not in VALID_BASES:
                raise ValueError(
                    f"Invalid base {base!r} at position {i} in plasmid "
                    f"{self.name!r}: only A, C, G, T are allowed"
                )
        if self.length == 0:
            raise ValueError(f"Plasmid {self.name!r} has an empty sequence")

    @classmethod
    def from_fasta(cls, path, circular=True):
        """Read one FASTA record from *path*.

        The header line starts with '>'. All following lines are
        concatenated, whitespace is stripped, and the sequence is
        upper-cased. Only the first record is read.
        """
        path = Path(path)
        name = None
        parts = []
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        # Only the first record is used.
                        break
                    # First word of the header is the name.
                    name = line[1:].strip().split()[0] if line[1:].strip() else path.stem
                else:
                    parts.append(line.upper())
        if name is None:
            raise ValueError(f"No FASTA record found in {path}")
        return cls(name, "".join(parts), circular=circular)

    def __repr__(self):
        return f"Plasmid(name={self.name!r}, length={self.length}, circular={self.circular})"
