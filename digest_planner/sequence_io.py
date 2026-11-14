"""Read plasmid sequences from FASTA files."""

# Only A, C, G, T are allowed in a sequence.
VALID_BASES = "ACGT"


class Plasmid:
    """One plasmid: a name and a DNA sequence.

    A plasmid is circular by default, so a recognition site may
    wrap around the end of the sequence back to the start.
    """

    def __init__(self, name, sequence, circular=True):
        self.name = name
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.circular = circular
        self.validate()

    def __str__(self):
        return self.name + " (" + str(self.length) + " bp)"

    def validate(self):
        """Reject any character other than A, C, G, T."""
        for i in range(len(self.sequence)):
            base = self.sequence[i]
            if base not in VALID_BASES:
                raise ValueError("Invalid base '" + base + "' at position "
                                 + str(i) + " in plasmid '" + self.name
                                 + "': only A, C, G, T are allowed")
        if self.length == 0:
            raise ValueError("Plasmid '" + self.name + "' is empty")

    @classmethod
    def from_fasta(cls, path, circular=True):
        """Read one FASTA record from a file.

        The first line starts with '>'. Every line after that is
        part of the sequence. Only the first record is read.
        """
        name = None
        parts = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line[0] == ">":
                    if name is not None:
                        break
                    words = line[1:].strip().split()
                    if words:
                        name = words[0]
                    else:
                        name = path.split("/")[-1]
                else:
                    parts.append(line.upper())
        if name is None:
            raise ValueError("No FASTA record found in " + path)
        return cls(name, "".join(parts), circular)
